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
import os
import time

data_path = '../data/'

# Read mesh
with io.XDMFFile(MPI.COMM_WORLD, data_path + 'mesh.xdmf', "r") as xdmf:
    domain = xdmf.read_mesh(name='Grid')

# Read meshtags (boundary information)
domain.topology.create_connectivity(domain.topology.dim-1, domain.topology.dim)
with io.XDMFFile(MPI.COMM_WORLD, data_path + 'mt.xdmf', "r") as xdmf:
    meshtags = xdmf.read_meshtags(domain, name="Grid")

# Define spaces
order_u = 1      # polynomial degree for displacement
order_p = 1      # polynomial degree for pressure
Ve = ufl.VectorElement("Lagrange", domain.ufl_cell(), order_u)
Pe = ufl.FiniteElement("Lagrange", domain.ufl_cell(), order_p)
state_space = fem.FunctionSpace(domain, Ve * Pe)

# Define functions and test functions
state = fem.Function(state_space)
test_state = ufl.TestFunction(state_space)

u, p = ufl.split(state)
v, q = ufl.split(test_state)

# Fibers
V_f = fem.VectorFunctionSpace(domain, ('CG', 1))       # space for fibers - continuous linear space
T_f = fem.TensorFunctionSpace(domain, ('CG', 1))       # space for outer product - continuous linear space

f0_array = np.loadtxt(data_path + 'fiber.txt')
s0_array = np.loadtxt(data_path + 'sheet.txt')
n0_array = np.loadtxt(data_path + 'sheetnormal.txt')

# pre compute outer product
outer_ff = np.zeros([*f0_array.shape, 3])
for i in range(len(f0_array)):
    outer_ff[i] = np.outer(f0_array[i], f0_array[i])
ff = fem.Function(T_f)
ff.vector.array = outer_ff.flatten()

f0 = fem.Function(V_f)    # In this tutorial fibers are defined in a linear topology
s0 = fem.Function(V_f)
n0 = fem.Function(V_f)

f0.vector.array = f0_array.flatten()
s0.vector.array = s0_array.flatten()
n0.vector.array = n0_array.flatten()

# Kinematics
dim = len(u)
I = ufl.variable(ufl.Identity(dim))
F = ufl.variable(I + ufl.grad(u))         # Deformation gradient
C = ufl.variable(F.T*F)                   # Right Green-Cauchy tensor

# Invariants
J = ufl.variable(ufl.det(F))
Ic = ufl.variable(ufl.tr(C))
Ic_bar = ufl.variable(J**(-2./3.) * Ic)
I4f = ufl.inner(C, ff)
I4s = ufl.inner(C*s0, s0)

# Define material law
a = fem.Constant(domain, PETSc.ScalarType(0.4))
b = fem.Constant(domain, PETSc.ScalarType(3.2))
a_f = fem.Constant(domain, PETSc.ScalarType(1.0))
b_f = fem.Constant(domain, PETSc.ScalarType(15.))
Psi_dev = a/(2.*b)*(ufl.exp(b*(Ic_bar-3.)) - 1.) \
    + a_f/(2.*b_f)*(ufl.exp(b_f*(I4f-1.)**2.) - 1.)

k = fem.Constant(domain, PETSc.ScalarType(1e2))
Psi_vol = p * (J - 1) - p**2 / (2*k)

# Get PK1 stress
P_dev = ufl.diff(Psi_dev,F)
P_vol = ufl.diff(Psi_vol,F)

# Define active stress
f = F*f0
lv_activation = fem.Constant(domain, PETSc.ScalarType(0.))
P_act = lv_activation * (ufl.outer(f,f0) + 0.2 * (F - ufl.outer(f,f0)))

# Boundary conditions
# Fixed base
V, _ = state_space.sub(0).collapse()
dofs = fem.locate_dofs_topological((state_space.sub(0), V), 2, meshtags.find(2))
u_fixed = fem.Function(V)
u_fixed.x.set(0.0)
bcs = [fem.dirichletbc(u_fixed, dofs, state_space.sub(0))]

# Endocardial pressure
normal_0 = ufl.FacetNormal(domain)
lv_pressure = fem.Constant(domain, PETSc.ScalarType(0.))
T = -lv_pressure * J * ufl.inv(F).T * normal_0

# Define integration domains
ds = ufl.Measure('ds', domain=domain, subdomain_data=meshtags)
dx = ufl.Measure("dx", domain=domain, metadata={'quadrature_degree': 4})

# Virtual work
deltaW_uu = ufl.inner(ufl.grad(v), P_dev)*dx + ufl.inner(ufl.grad(v), P_act)*dx - ufl.inner(v, T)*ds(3)
deltaW_up = ufl.inner(ufl.grad(v), P_vol)*dx
deltaW_pu = q * (J - 1.) * dx
deltaW_pp = - 1/k * p * q * dx
deltaW = deltaW_uu + deltaW_up + deltaW_pu + deltaW_pp


# Setting up the problem and solver
problem = fem.petsc.NonlinearProblem(deltaW, state, bcs)
solver = nls.petsc.NewtonSolver(domain.comm, problem)
solver.atol = 1e-8
solver.rtol = 1e-8
solver.error_on_nonconvergence = False

# Read data for activation and pressure
activation = np.loadtxt(data_path + 'activation.txt')
pressure = np.loadtxt(data_path + 'pressure.txt')

from scipy.interpolate import interp1d
func_act = interp1d(activation[:,0], activation[:,1])
func_pres = interp1d(pressure[:,0], pressure[:,1])


# First step: inflate LV to initial pressure
# Initializing files to save results (clean first)
if os.path.isfile('inflation.xdmf'):
    os.remove('inflation.xdmf')
    os.remove('inflation.h5')

disp_file = io.XDMFFile(MPI.COMM_WORLD, 'inflation.xdmf', 'w')
disp_file.write_mesh(domain)

# Incremmentally load the lv using 10 timesteps
nsteps = 21
pres_values = np.linspace(0, func_pres(0), nsteps+1)[1:]
for i in range(len(pres_values)):
    lv_pressure.value = pres_values[i]
    num_its, converged = solver.solve(state)
    if not converged:
        disp_file.close()
        assert converged
    u, p = state.split()
    disp_file.write_function(u, i)
    print(f"Time step {i}, Number of iterations {num_its}, pressure {lv_pressure.value}")
disp_file.close()

inf_state = fem.Function(state_space)
inf_state.vector.array = np.copy(state.vector.array)

# Second step: Cardiac cycle
state.vector.array = np.copy(inf_state.vector.array)

if os.path.isfile('cycle.xdmf'):
    os.remove('cycle.xdmf')
    os.remove('cycle.h5')
disp_file = io.XDMFFile(MPI.COMM_WORLD, 'cycle.xdmf', 'w')
disp_file.write_mesh(domain)

start = time.time()
nsteps = 501
times = np.linspace(0, 1, nsteps)[1:]    # We already solved the first step
for i, t in enumerate(times):
    # Set pressure/activation values
    lv_pressure.value = func_pres(t)
    lv_activation.value = func_act(t)

    # Solve
    num_its, converged = solver.solve(u)

    if not converged:
        disp_file.close()
        assert converged

    # Save output
    u, p = state.split()
    disp_file.write_function(u, i)

    print(f"Time step {i}, Number of iterations {num_its}, pressure {lv_pressure.value}, activation {lv_activation.value}")
disp_file.close()

end = time.time()

print(f"Simulation time {end-start} s")

f0.vector.destroy()
s0.vector.destroy()
n0.vector.destroy()
u.vector.destroy()
