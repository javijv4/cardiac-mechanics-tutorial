# cardiac-mechanics-tutorial
A cardiac mechanics tutorial using fenics-x

## Getting Started
Before starting, install conda or miniconda (either is fine, miniconda is a lightweight version of conda) following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once conda is installed, create an environment by running the following,
```
conda create -n cardiac-mechanics python=3.10
```
Here, `cardiac-mechanics` is the name of the environment. Now activate the environment and install fenics-x and all the dependencies required to run the tutorial (basically, fenics-x and jupyter). This step might take a while, depending on your system.
```
conda activate cardiac-mechanics
conda install -c conda-forge fenics-dolfinx mpich jupyter
conda install scipy matplotlib
```
Clone the repository (this will make a folder from wherever folder you are calling the command):
```
git clone https://github.com/javijv4/cardiac-mechanics-tutorial
```
To run the second tutorial, an specialized package [ambit](https://github.com/marchirschvogel/ambit) is needed. This can be installed using (make sure the cardiac-mechanics environment is activated),
```
python -m pip install ambit-fe==1.0.7
```

## Tutorials
Before running any tutorial make sure to activate the environment,
```
conda activate cardiac-mechanics
```
Each tutorial has a `python` file (`.py`) and a `jupyter-notebook` version (`.ipynb`). The code is the same, but the `jupyter-notebook` has more information about the steps.

### Tutorial 1 - Simple Left Ventricle
This tutorial consists in an idealized left ventricle geometry. For the tutorial purposes a very coarse discretization is used. Chamber pressure and cardiac muscle activation curves are given to recreate the cardiac cycle. The provided `python` and `jupyter-notebook` files detail the setup of the problem, from reading the input files, the hyperelastic formulation, boundary conditions and solver setup.

### Tutorial 2 - Coupled Left Ventricle
In this tutorial the same setup as Tutorial 1 is used but instead of prescribing a chamber pressure, the solid 3D model is coupled with a 0D representation of the cardiovascular system. This is a more complex setup, so instead of directly using fenics-x to set up the problem, we use the package [ambit](https://github.com/marchirschvogel/ambit) that is also based on fenics-x. The provided `python` and `jupyter-notebook` files detail the setup of the file and options required to setup the 3D-0D problem.

## Resources
For more examples of how to use fenics-x, we recommend looking at this comprehensive [tutorial](https://jsdokken.com/dolfinx-tutorial/).

The cardiac mechanics implementation in Tutorial 1 is loosely based on the following repositories: [fenicsx-pulse](https://github.com/finsberg/fenicsx-pulse) and [ambit](https://github.com/marchirschvogel/ambit).
