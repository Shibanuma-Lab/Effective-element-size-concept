C==============================================================================C
C                 Umat for ductile fracture modeling with                      C
C                  Stress Modified Fracture strain model                       C
C          Developed by Shibanuma Lab from The University of Tokyo             C
C==============================================================================C
C                                                                              C
C                                                                              C
C                                                                              C
C **************** Hardening_Index - Indicator of stage ************************
C      0: Elastic stage                                                        *
C      1: Hardening stage                                                      *
C      2: Softening stage                                                      *
C      3: Fracture stage                                                       *
C -----------------Frac_Index- Indicator of stage------------------------------*
C      0: Other stage                                                          *
C      1: Fracture stage                                                       *
C------------------   Input parameters   --------------------------------------*
C PROPS(1)   - Young's modulus                                                 *
C PROPS(2)   - Poisson's ratio                                                 *
C PROPS(3)   - parameter A                                                     *
C PROPS(4)   - parameter B                                                     *
C PROPS(5..) - Mises hardening data                                            *
C                                                                              *
C NHARD  - Number of hardening data                                            *
C SLOPE  - Slope of softening curve                                            *
C SRMIN  - Ratio of minimum stress to maximum stress                           *
C NIPD    - Number of element delete judge                                     *
C Damage_Index - damage indicator (0:damaged, 1: active)                       *
C ******************************************************************************
C


      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION TIME(2)
#include <SMAAspUserSubroutines.hdr>
	  INTEGER Gausspoint_Act_Index
	  pointer(ptr_Gausspoint_Act_Index,Gausspoint_Act_Index(200000,8))
	  IF(lop.EQ.0) THEN
		 ptr_Gausspoint_Act_Index = SMAIntArrayCreate(1,200000*8,1)
	  END IF
      RETURN
      END SUBROUTINE UEXTERNALDB


      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'
#include <SMAAspUserSubroutines.hdr>
      CHARACTER*8 CMNAME
	  INTEGER Gausspoint_Act_Index
	  pointer(ptr_Gausspoint_Act_Index,Gausspoint_Act_Index(200000,8))

      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 PREDEF(1),DPRED(1),PROPS(NPROPS), COORDS(3),DROT(3, 3),
     3 DFGRD0(3,3),DFGRD1(3,3),EELAS(6),EPLAS(6),FLOW(6),HARD(3),
     4 SOFT(2,2000),DEVSTRESS(6),TIME(2)
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

      NHARD=5346
      SLOPE=-5000.0D0
      SRMIN=0.02D0
      NIPD=4
C
C---------------------------------------------------------------------
C Variables from previous increment
      CALL ROTSIG(STATEV(1),DROT,EELAS,2,NDI,NSHR)          ! 读取弹性应变分量，塑性应变分量，
      CALL ROTSIG(STATEV(NTENS+1),DROT,EPLAS,2,NDI,NSHR)
	  ptr_Gausspoint_Act_Index = SMAIntArrayAccess(1)
	   
	  EQPLAS=STATEV(1+2*NTENS)                              ! 读取等效塑性应变
      DAMAGE=STATEV(2+2*NTENS)
      SYIELD=STATEV(3+2*NTENS)
      FFLAG=1
      SYIELDMAX=STATEV(3*NTENS)
      Hardening_Index=STATEV(3+3*NTENS)
	  SEQPLAS=STATEV(5+3*NTENS)
      Damage_Index=STATEV(4*NTENS)
      Frac_Index=STATEV(1+4*NTENS)    
      JSLOPE=STATEV(2+4*NTENS)
      
C  Elastic properties
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2), ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)                             !   3K , 3k = E /(1-2ν )
      EG2=EMOD/(ONE+ENU)                                    !   2G , 2G = E /(1+υ)
      EG=EG2/TWO                                            !   G  , G = E /2(1+υ)
      EG3=THREE*EG                                          !   3G
      ELAM=(EBULK3-EG2)/THREE                               !   λ  ,  λ = (3k-2G)/3

C Elastic stiffness (elastic D matrix)	 
      DO K1=1, NDI
         DO K2=1, NDI
            DDSDDE(K2, K1)=ELAM
         END DO
         DDSDDE(K1, K1)=EG2+ELAM
      END DO
      
      DO K1=NDI+1, NTENS
         DDSDDE(K1, K1)=EG
      END DO

C Calculation of trial stress and elastic strain tensors    !按完全弹性变形计算试算应力Δσ
      DO K1=1, NTENS
         DO K2=1, NTENS
            STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
         END DO
         EELAS(K1)=EELAS(K1)+DSTRAN(K1)                     ! 弹性应变分量
      END DO
C Equivalent Von Mises stress
	  SMISES=(STRESS(1)-STRESS(2))**TWO+(STRESS(2)-STRESS(3))**TWO
     1 +(STRESS(3)-STRESS(1))**TWO
      DO K1=NDI+1,NTENS
         SMISES=SMISES+SIX*STRESS(K1)**TWO
      END DO
      SMISES=SQRT(SMISES/TWO)    

C Flow directions and Mean stress	                        !计算流动方向  
	  SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
      DO K1=1,NDI
         FLOW(K1)=(STRESS(K1)-SHYDRO)/SMISES
      END DO
      DO K1=NDI+1, NTENS
         FLOW(K1)=STRESS(K1)/SMISES
      END DO
	  
C---------------------------------------------------------------------
C=============================Start===================================
C---------------------------------------------------------------------
C
      Hardening_Index=0  
      Damage_Index=1 
      IF(DAMAGE.LT.ONE) THEN
C---------------------------------------------------------------------
C 0 - Elastic stage
C---------------------------------------------------------------------    
      Hardening_Index=0
      Damage_Index=1   
      CALL Y_SUB(EQPLAS,SYIELD,HARD,PROPS,NPROPS,NHARD)        
C---------------------------------------------------------------------
C 1 - Hardening stage
C---------------------------------------------------------------------   
 	    IF(SMISES.GT.(ONE+TOLER)*SYIELD) THEN
            Hardening_Index=1
	        DEQPL=ZERO
            DO KEWTON=1,NEWTON
               RHS=SMISES-EG3*DEQPL-SYIELD
               DEQPL=DEQPL+RHS/(EG3+HARD(1))
               SYIELD=ZERO
               HARD(1)=ZERO    !????没用吧  
               CALL Y_SUB(EQPLAS+DEQPL,SYIELD,HARD,PROPS,NPROPS,NHARD)      
               IF(ABS(RHS).LE.TOLER*SYIELD) GOTO 12
            END DO
   12 CONTINUE
        ENDIF	       
      ELSE
C---------------------------------------------------------------------
C 2 - Softening stage
C---------------------------------------------------------------------    
        Damage_Index=0      
	    IF(Frac_Index.NE.1) THEN
            IF(JSLOPE.EQ.0) THEN
                SYIELDMAX=SYIELD
                SEQPLAS=EQPLAS
                JSLOPE=1
            ENDIF      
C Softening curve with transition	
            SRTRN=0.019D0	  
            NSOFT=100	  
            IF(SRTRN.GT.SRMIN) THEN
                NS=NSOFT
                E1=SEQPLAS
                S1=SYIELDMAX
                S2=SRTRN*SYIELDMAX
                S3=SRMIN*SYIELDMAX
                E2=E1+(S2-S1)/SLOPE
	            
                SOFT(1,1)=0.99D0*E1
                SOFT(2,1)=S1-0.01D0*SLOPE
                SOFT(1,2)=E1
                SOFT(2,2)=S1
                SOFT(1,3)=E2
                SOFT(2,3)=S2
	            
                E3=E2-TWO*(S2-S3)/SLOPE
                B=(SLOPE/TWO)/(E2-E3)
                DO K1=4,NS-1
                   SOFT(1,K1)=E2+(K1-2)/(NS-3)*(E3-E2)
                   SOFT(2,K1)=S3+B*(SOFT(1,K1)-E3)**TWO
                ENDDO
	            
                SOFT(1,NS)=100
                SOFT(2,NS)=S3
            ELSE
C Softening curve without transition (only linear part)
                NS=4
                E1=SEQPLAS
                S1=SYIELDMAX
                S3=SRMIN*SYIELDMAX
                E3=E1+(S3-S1)/SLOPE	  
			    
                SOFT(1,1)=0.99D0*E1
                SOFT(2,1)=S1-0.01D0*SLOPE
                SOFT(1,2)=E1
                SOFT(2,2)=S1
                SOFT(1,3)=E3
                SOFT(2,3)=S3
                SOFT(1,4)=100
                SOFT(2,4)=S3     
            ENDIF

									  
						 
						 
						 
						  
            CALL S_SUB(EQPLAS,SYIELD,HARD,SOFT,NS)	  
	        IF(SMISES.GT.(ONE+TOLER)*SYIELD) THEN
	            Hardening_Index=2
	            DEQPL=ZERO
	            DO KEWTON=1,NEWTON
	               RHS=SMISES-EG3*DEQPL-SYIELD
	               DEQPL=DEQPL+RHS/(EG3+HARD(1))	               
	               SYIELD=ZERO
                   HARD(1)=ZERO
			       CALL S_SUB(EQPLAS+DEQPL,SYIELD,HARD,SOFT,NS)
			       IF(ABS(RHS).LE.TOLER*SYIELD) GOTO 15
                END DO
   15 CONTINUE
			    IF(SYIELD.LE.SRMIN*SYIELDMAX*1.001D0) THEN
                   Frac_Index=1
                   Hardening_Index=3
                ENDIF      
            ENDIF
        ELSE              
C---------------------------------------------------------------------
C 3 - Fractured stage
C---------------------------------------------------------------------
			 SYIELD=SRMIN*SYIELDMAX
             IF(SMISES.GT.(ONE+TOLER)*SYIELD) THEN
			    Hardening_Index=3
                DEQPL=ZERO
                HARD(1)=ZERO
			    DO KEWTON=1,NEWTON
                   RHS=SMISES-EG3*DEQPL-SYIELD
                   DEQPL=DEQPL+RHS/EG3
                   IF(ABS(RHS).LE.TOLER*SYIELD) GOTO 25
                END DO
   25 CONTINUE
             ENDIF	      
         ENDIF
      ENDIF	      
      
C---------------------------------------------------------------------
C=================Update of stress, strain and damage=================
C---------------------------------------------------------------------
      IF(Hardening_Index.GE.1) THEN
        DO K1=1,NDI
           STRESS(K1)=FLOW(K1)*SYIELD+SHYDRO
           EPLAS(K1)=EPLAS(K1)+THREE/TWO*FLOW(K1)*DEQPL
           EELAS(K1)=EELAS(K1)-THREE/TWO*FLOW(K1)*DEQPL
        END DO        
        DO K1=NDI+1,NTENS
           STRESS(K1)=FLOW(K1)*SYIELD
           EPLAS(K1)=EPLAS(K1)+THREE*FLOW(K1)*DEQPL
           EELAS(K1)=EELAS(K1)-THREE*FLOW(K1)*DEQPL
        END DO      
        EQPLAS=EQPLAS+DEQPL
        SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
        DO K1=1, NDI
           DEVSTRESS(K1)=STRESS(K1)-SHYDRO
        END DO
        DO K1=NDI+1,NTENS
           DEVSTRESS(K1)=STRESS(K1)
        END DO 
C---------------------------------------------------------------------
C Damage based on trixiality (and Lode angle parameter)
        TRIX=SHYDRO/SYIELD
	    J2=((DEVSTRESS(1)-DEVSTRESS(2))**TWO+(DEVSTRESS(1)-DEVSTRESS(3))**TWO+(DEVSTRESS(3)-DEVSTRESS(2))**TWO)/6
     1   +DEVSTRESS(4)**TWO+DEVSTRESS(5)**TWO+DEVSTRESS(6)**TWO	  
        J3=DEVSTRESS(1)*DEVSTRESS(2)*DEVSTRESS(3)+TWO*DEVSTRESS(4)*DEVSTRESS(5)
     1   *DEVSTRESS(6)-DEVSTRESS(1)*(DEVSTRESS(6))**TWO-DEVSTRESS(2)*
     2   (DEVSTRESS(5))**TWO-DEVSTRESS(3)*(DEVSTRESS(4))**TWO

	    ALODEpara=THREE*(THREE**(ONE/TWO))*J3/(TWO*J2**(THREE/TWO))
        IF (ALODEpara.LE.-ONE ) THEN
	       ALODEpara=-ONE
        ENDIF	  
        IF (ALODEpara.GE.ONE ) THEN
	       ALODEpara=ONE
        ENDIF	  	  
	    ALODE=(ONE/THREE)*ACOS(ALODEpara)
	    ALODEnorma=ONE-TWO*THREE*ALODE/ACOS(-ONE)

		MINTRA= ONE/THREE
        IF (TRIX.LE.MINTRA) THEN
	       CRIEQPLAS=PROPS(3)*EXP(-PROPS(4)*TRIX)
        ELSE 
	       CRIEQPLAS=PROPS(3)*EXP(-PROPS(4)*TRIX)
        ENDIF		 	  	 
        DDAMAGE=DEQPL/CRIEQPLAS
        DAMAGE=DAMAGE+DDAMAGE  
        IF(DAMAGE.GE.ONE) THEN
          DAMAGE=ONE
          Damage_Index=0
        ENDIF  
C---------------------------------------------------------------------
C Jacobian (Material tangent)
        EFFG=EG*SYIELD/SMISES
        EFFG2=TWO*EFFG
        EFFG3=THREE/TWO*EFFG2
        EFFLAM=(EBULK3-EFFG2)/THREE
        EFFHRD=EG3*HARD(1)/(EG3+HARD(1))-EFFG3
        DO K1=1, NDI
            DO K2=1, NDI
                DDSDDE(K2, K1)=EFFLAM
            END DO
            DDSDDE(K1, K1)=EFFG2+EFFLAM
        END DO
        DO K1=NDI+1, NTENS
            DDSDDE(K1, K1)=EFFG
        END DO	    
        DO K1=1, NTENS
            DO K2=1, NTENS
                DDSDDE(K2, K1)=DDSDDE(K2, K1)+EFFHRD*FLOW(K2)*FLOW(K1)
            END DO
         END DO
      ENDIF 
C---------------------------------------------------------------------
C State variables       
           
      DO K1=1,8
	     IF (NPT.EQ.K1) Gausspoint_Act_Index(NOEL,K1)=Damage_Index
	     Damage_IndexSUM=Damage_IndexSUM+Gausspoint_Act_Index(NOEL,K1)
      END DO
	  IF (Damage_IndexSUM.LE.8-NIPD) THEN
	     FFLAG=0
	  ENDIF
      
      DO K1=1, NTENS
		 STATEV(K1)=EELAS(K1)
		 STATEV(K1+NTENS)=EPLAS(K1)
      END DO
      STATEV(1+2*NTENS)=EQPLAS
      STATEV(2+2*NTENS)=DAMAGE
      STATEV(3+2*NTENS)=SYIELD
      STATEV(4+2*NTENS)=FFLAG
      STATEV(5+2*NTENS)=DEQPL
      STATEV(3*NTENS)=SYIELDMAX
      STATEV(1+3*NTENS)=TRIX
      STATEV(2+3*NTENS)=ALODE
      STATEV(3+3*NTENS)=Hardening_Index
      STATEV(4+3*NTENS)=CRIEQPLAS
      STATEV(5+3*NTENS)=SEQPLAS
      STATEV(4*NTENS)=Damage_Index
      STATEV(1+4*NTENS)=Frac_Index
      STATEV(2+4*NTENS)=JSLOPE

      STATEV(3+4*NTENS)=ALODEpara
      STATEV(4+4*NTENS)=ALODEnorma
      RETURN
      END
C
C **************************************************************************** C
C        Y_SUB and S_SUB are used for hardening and yield stress update        C
C        S_SUB Used For Linear Softening                                       C
C **************************************************************************** C
C
      SUBROUTINE Y_SUB(EQPLAS,SYIEL0,HARD,PROPS,NPROPS,NHARD)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION PROPS(NPROPS), HARD(3)
C
      IF(EQPLAS.GE.PROPS(NPROPS)) THEN
        HARD(1) = ZERO
        SYIEL0  = PROPS(NPROPS-1)		
      ELSE
		  DO K1=1,NHARD
		  EQPL1=PROPS(2*K1+6)
		  IF(EQPLAS.LT.EQPL1) THEN
			EQPL0=PROPS(2*K1+4)
			SYIEL00=PROPS(2*K1+3)
			SYIEL10=PROPS(2*K1+5)
			HARD(1)=(SYIEL10-SYIEL00)/(EQPL1-EQPL0)
			SYIEL0=SYIEL00+(EQPLAS-EQPL0)*HARD(1)
			GOTO 91
		  ENDIF
		  END DO
	  ENDIF
   91 CONTINUE
      
      RETURN
      END
C
      SUBROUTINE S_SUB(EQPLAS,SYIEL0,HARD,SOFT,NS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SOFT(2,NS), HARD(3)
C
      IF(EQPLAS.GE.SOFT(1,NS)) THEN
        HARD(1) = ZERO
        SYIEL0  = SOFT(2,NS)
      ELSE
        DO K1=1,NS-1
        EQPL1=SOFT(1,K1+1)
        IF(EQPLAS.LT.EQPL1) THEN
          EQPL0=SOFT(1,K1)
          SYIEL00=SOFT(2,K1)
          SYIEL10=SOFT(2,K1+1)
		  HARD(1)=(SYIEL10-SYIEL00)/(EQPL1-EQPL0)
          SYIEL0=SYIEL00+(EQPLAS-EQPL0)*HARD(1)
          GOTO 92
        ENDIF
        END DO
      ENDIF
   92 CONTINUE
      
      RETURN
      END



