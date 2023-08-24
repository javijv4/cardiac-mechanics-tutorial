# cardiac-mechanics-tutorial
A cardiac mechanics tutorial using fenicsx

## Getting Started
Before start, we need to install conda or miniconda (either is fine, miniconda is a lightweight version of conda) following the instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once conda is installed, we first need to create an environment by running the following,
```
conda create -n cardiac-mechanics
```
here, `cardiac-mechanics` is the name of the environment. Now we need to install fenics-x and all the dependencies required to run the tutorial (basically, fenics-x and jupyter). This step might take a while depending on your system. 
```
conda install -c conda-forge fenics-dolfinx h5py pyvista -y
conda install -c conda-forge jupyter
```
Now that all the packages are installed, we need to activate the environment,
```
cd cardiac-mechanics-tutorial
```
And we proceed with cloning the repository (this will make a folder from wherever folder you are calling the command):
```
git clone https://github.com/javijv4/cardiac-mechanics-tutorial
```
Enter the newly created folder and run the jupyter-notebook.
```
cd cardiac-mechanics-tutorial
jupyter-notebook tutorial.ipynb
```
