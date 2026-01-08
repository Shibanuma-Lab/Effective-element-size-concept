## Automated Effective Element Size–dependent Parameter Identification for Ductile Fracture Modeling (Abaqus–Python Framework)

This repository provides a fully automated framework for effective element size–dependent parameter identification of a ductile fracture model developed by Shibanuma Lab, The University of Tokyo.  
It implements a hybrid experimental–numerical strategy, further extended by introducing the **effective element size (s)** as a control variable, integrating Python scripting with Abaqus analyses to automate the entire identification process — including input generation, job execution, data extraction, and error evaluation.  
Through this framework, key model parameters are systematically identified, including:
- hardening exponent `n`
- fracture strain parameter `fs` (identified at multiple element sizes)
- material parameters `A` and `B` (identified at multiple element sizes)
- material parameters `A(s)` and `B(s)`, expressed as functions of the effective element size


## Framework Overview
The framework consists of two main Python modules:

- **ALLEXE.py** — Serves as the execution script for running the entire identification procedure.  
  This file defines user-specific settings such as 
  - the reference element size and additional target element sizes,
  - parameter search ranges,
  - input file names,
  - computational resources and parallel job options.
  - other customizable options
  
- **ALLidentifcation.py** — Contains the main functions and core procedures for automated parameter identification.  
  This file **normally does not require modification**. It defines the underlying logic of the identification process, including routines for input preparation, simulation control, and error evaluation.


## Methodology
The entire process is organized into the following sequential stages:
1. **Element Size Configuration (preparing stage)** 
   Prepare mesh templates at reference element size and a set of target sizes {`s₀₁`,..., `s₀ᵢ`} for both tensile and Arcan tests.
2. **Hardening identification (n stage)**  
   The hardening exponent `n` is identified once at the reference element size, which remains constant throughout the entire process.
3. **Fracture strain identification (fs stage)** 
   The fracture strain parameter `fs` is identified at each target element size with fixed `n`.
4. **Material parameter identification (A–B stage)** 
   Material parameters `A` and `B` are identified at each target element size with fixed `n` and `fs(s₀ᵢ)`.
5. **Evaluation of effective element size (s stage)**  
   For each target element size, an effective element size `sᵢ` is evaluated based on the identified parameters.
6. **Effective element size-dependent parameter construction**  
  The material parameters `A` and `B` are represented as functions of the effective element size.

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
  (e.g., `Tensile_s0.inp`, `Arcan_s0.inp`, `Tensile_s1.inp`, `Arcan_s1.inp`).

- **User subroutine file** — Fortran file implementing the constitutive model to be compiled by Abaqus  
  (e.g., `UMAT_name.f`).

- **Batch file for Abaqus execution** — Template for running Abaqus jobs via the command line  
  (e.g., `make023.bat`).

All these files must be placed in the same directory as the Python scripts.


### 3. Running the Framework

All parameter and path settings are defined in **`ALLEXE.py`**, which serves as the main entry point of the framework.  
In this file, users can specify:
- Reference and target element sizes,
- Material name and input file names  
- Parameter search ranges for `n`, `fs`, and `B`  
- Number of CPUs and parallel job settings  

Once configured, execute the following command in a terminal:
```bash
python ALLEXE.py
```

### 4. Output Files

The outputs record parameter evolution, error values, and identification results at different element sizes, as well as the final functions.
For parameters identified at multiple element sizes, the results are grouped into subfolders corresponding to each element size.

- **ALLResults/**  
  Contains all final result files, including the summary `.csv` files. These files record the optimization history and the final identified values.

- **n/**, **fs/**, **b/**  
  Each folder stores intermediate data, temporary files, and Abaqus results corresponding to the three stages of identification.
  
- **ALLinput/**  
  Includes all input-related data such as experimental datasets.

In addition, a dedicated `functions.csv` file is generated to store the functions of the material parameters `A` and `B` based on the effective element size.


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

