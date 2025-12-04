# Automated Parameter Identification for Ductile Fracture Modeling (Abaqus–Python Framework)

This project provides a fully automated framework for parameter identification of a ductile fracture model developed by Shibanuma Lab, The University of Tokyo.  
It implements a hybrid experimental–numerical strategy, integrating Python scripting with Abaqus analyses to automate the entire identification process — including input generation, job execution, data extraction, and error evaluation.  
Within this framework, Abaqus simulations are incorporated into a numerical optimization loop that minimizes the discrepancy between experimental and simulated responses.

Through this framework, key model parameters are systematically identified, including:
- Hardening exponent `n`
- Fracture strain `fs`
- Material constants `A` and `B`


## Framework Overview
The framework consists of two main Python modules:

- **ALLEXE.py** — Serves as the execution script for running the entire identification procedure.  
  This file defines user-specific settings such as parameter ranges, input file names, and other customizable options required for each identification stage.
- **ALLidentifcation.py** — Contains the main functions and core procedures for automated parameter identification.  
  This file **normally does not require modification**. It defines the underlying logic of the identification process, including routines for input preparation, simulation control, and error evaluation.

The entire process is organized in three sequential stages:
1. **Hardening identification (n stage)**  
2. **Fracture strain identification (fs stage)** 
3. **Material parameter identification (A–B stage)** 

All results and simulation data are automatically saved in `.csv` and `.pkl` files for restart and post-processing.


## Usage

### 1. Environment Requirements

- **Abaqus** with user subroutine (UMAT) capability  
- **Python environment** with libraries installed:  
  `numpy`, `pandas`, `scipy`, `matplotlib`, `sympy`, `mpmath`  
- Ensure that Abaqus supports **batch execution of user subroutine (UMAT) jobs** through the command line.


### 2. Input Preparation

Before running the framework, prepare the following input data and files in the working directory:

- **Experimental data files** — Contain tensile and Arcan test results, including DIC-measured data  
  (e.g., `material_ss.csv`, `Arcan_Load.csv`).

- **Abaqus input templates** — Finite element models for tensile and Arcan tests  
  (e.g., `Tensile.inp`, `Arcan.inp`).

- **User subroutine file** — Fortran file implementing the constitutive model to be compiled by Abaqus  
  (e.g., `UMAT_name.f`).

- **Batch file for Abaqus execution** — Template for running Abaqus jobs via the command line  
  (e.g., `make023.bat`).

All these files must be placed in the same directory as the Python scripts.


### 3. Running the Framework

All parameter and path settings are defined in **`ALLEXE.py`**, which serves as the main entry point of the framework.  
In this file, users can specify:
- Material name and input file names  
- Parameter search ranges for `n`, `fs`, and `B`  
- Number of CPUs and parallel job settings  

Once configured, execute the following command in a terminal:
```bash
python ALLEXE.py
```

### 4. Output Files

After the identification process is completed, the framework automatically organized all result and intermediate files into subfolders.  
These files record parameter evolution, error values, and the best-fit results for each stage.

- **ALLResults/**  
  Contains all final result files, including the summary `.csv` files. These files record the optimization history and the final identified values.

- **n/**, **fs/**, **b/**  
  Each folder stores intermediate data, temporary files, and Abaqus results corresponding to the three stages of identification.
  
- **ALLinput/**  
  Includes all input-related data such as experimental datasets.
  
All results are consolidated in the `Results_ALL.csv` files located in the working directory.


### 5. Notes

- **Path configuration**  
  All files are read and written relative to the working directory where the scripts are executed.  

- **Abaqus version compatibility**  
  Commands in the batch file (e.g., `make023.bat`, `ab023.bat`) should be consistent with your Abaqus installation name  
  (for example, replace `ab2023` with `abaqus2021` if necessary).

- **Restart mechanism**  
  If a previous run was interrupted, the framework automatically detects existing result files and continue without recalculating completed jobs.
  Please ensure that no corrupted `.dat` files are present in the working directory before restarting.

- **Mesh type**  
  The Arcan model in this framework uses a **2D (plane) mesh** by default. If a 3D mesh is used instead, the call to the `ArcanMesh` function in the Python script can be commented out.  


