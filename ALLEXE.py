import os
from ALLidentifcation import ALL_START_NFsB
from ALLidentifcation import initialize_parameters
from ALLidentifcation import initialize_parameters_all

current_directory = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_directory)
# Default experimental settings; can be modified manually
thickness, nz, dstep1= initialize_parameters_all() #
ftstep2n, dtstep2n, dstep2n= initialize_parameters("n")
ftstep2fs, dtstep2fs, dstep2fs= initialize_parameters("fs")
ftstep2b, dtstep2b, dstep2b= initialize_parameters("b")

# Set base parameters
material = "material"
inpname_tensile = "Tensile"
inpname_arcan = "Arcan"
UMATname = "UMAT_name"
# Set the range for the hardening exponent N
ln = 0.0
un = 0.2
# Set the range for the fracture strain fs
lfs = 0.0
ufs = 1.5
# Set the range for the material dependant parameters B
lb = 0.
ub = 4
# Set additional parameters
cpus = 2
jobnumber = 3  # Select 1 or 3; controls the number of parallel Abaqus runs
# ---------------------Starting Calculation---------------------
ALL_START_NFsB(material, inpname_tensile, inpname_arcan, UMATname, dstep1,
                    ftstep2n, dtstep2n, dstep2n, un, ln,
                    ftstep2fs, dtstep2fs, dstep2fs, ufs, lfs,
                    ftstep2b, dtstep2b, dstep2b, ub, lb,
                    nz, thickness, cpus, jobnumber)
