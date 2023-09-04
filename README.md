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

## Running the tutorial
Every time you want to run the tutorial after the first time, you need to activate the environment and call jupyter.
```
conda activate cardiac-mechanics
jupyter-notebook tutorial.ipynb
```
