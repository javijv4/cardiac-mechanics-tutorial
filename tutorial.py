#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 17:05:37 2023

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import io, fem, nls

data_path = 'data/'

def visualize_function(fname, u):
    mesh = u.function_space.mesh
    with io.XDMFFile(mesh.comm, fname, "w") as file:
        file.write_mesh(mesh)
        file.write_function(u)


# Read mesh
with io.XDMFFile(MPI.COMM_WORLD, data_path + 'mesh.xdmf', "r") as xdmf:
    mesh = xdmf.read_mesh(name='Grid')

# Read meshtags (boundary information)
mesh.topology.create_connectivity(mesh.topology.dim-1, mesh.topology.dim)
with io.XDMFFile(MPI.COMM_WORLD, data_path + 'mt.xdmf', "r") as xdmf:
    meshtags = xdmf.read_meshtags(mesh, name="Grid")


# Defining space functions
V_u = fem.VectorFunctionSpace(mesh, ('CG', 1))       # space for displacement - continuous quadratic space
V_f = fem.VectorFunctionSpace(mesh, ('CG', 1))       # space for fibers - continuous linear space

# Define test and trial function
du = ufl.TrialFunction(V_u)            # Incremental displacement
var_u = ufl.TestFunction(V_u)             # Test function
u = fem.Function(V_u, name="Displacement")  # Function to store solution


# Read fibers
f = fem.Function(V_f)    # In this tutorial fibers are defined in a linear topology
s = fem.Function(V_f)
n = fem.Function(V_f)

f.vector.array = np.loadtxt(data_path + 'fiber.txt').flatten()
s.vector.array = np.loadtxt(data_path + 'sheet.txt').flatten()
n.vector.array = np.loadtxt(data_path + 'sheetnormal.txt').flatten()

# Define surface normal direction
n0 = ufl.FacetNormal(mesh)

# Kinematics
dim = len(u)
I = ufl.variable(ufl.Identity(dim))
F = ufl.variable(I + ufl.grad(u))         # Deformation gradient
C = ufl.variable(F.T*F)                   # Right Green-Cauchy tensor

# Invariants
Ic = ufl.variable(ufl.tr(C))
IIIc = ufl.variable(ufl.det(C))
Ic_bar = ufl.variable(IIIc**(-1./3.) * Ic)

# Define material law
mu = fem.Constant(mesh, PETSc.ScalarType(10))
Psi_dev = (mu/2.) * (Ic_bar - 3.)
K = fem.Constant(mesh, PETSc.ScalarType(100))
Psi_vol = (K/2.) * (IIIc - 1.)**2

# Get PK1 stress
P = 2.*ufl.diff(Psi_dev,F) + 2.*ufl.diff(Psi_vol,F)

# Boundary conditions
# Fixed base
u_bc = np.array((0,) * mesh.geometry.dim, dtype=PETSc.ScalarType)
base_dofs = fem.locate_dofs_topological(V_u, meshtags.dim, meshtags.find(2))
bcs = [fem.dirichletbc(u_bc, base_dofs, V_u)]

lv_pressure = fem.Constant(mesh, PETSc.ScalarType(1))
T = -lv_pressure*n0

# Define surfaces
ds = ufl.Measure('ds', domain=mesh, subdomain_data=meshtags)
dx = ufl.Measure("dx", domain=mesh)

# Define residual
F = ufl.inner(ufl.grad(var_u), P)*dx
F += -ufl.inner(var_u, T)*ds(3)

problem = fem.petsc.NonlinearProblem(F, u, bcs)

solver = nls.petsc.NewtonSolver(mesh.comm, problem)

# Set Newton solver options
solver.atol = 1e-8
solver.rtol = 1e-8
solver.convergence_criterion = "incremental"

num_its, converged = solver.solve(u)

f.vector.destroy()
s.vector.destroy()
n.vector.destroy()
u.vector.destroy()