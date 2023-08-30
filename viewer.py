#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 16:46:32 2023

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
from dolfinx import plot, fem
import pyvista
pyvista.set_jupyter_backend('trame')
pyvista.global_theme.trame.server_proxy_enabled
pyvista.global_theme.trame.server_proxy_prefix

plot_in_line = True

def plot_mesh(mesh):
    if not plot_in_line: return
    tdim = mesh.topology.dim
    pyvista.start_xvfb()
    p = pyvista.Plotter()
    topology, cell_types, geometry = plot.create_vtk_mesh(mesh, tdim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    p.add_mesh(grid, show_edges=True)
    p.show()


def plot_meshtags(mesh, meshtags):
    if not plot_in_line: return
    tdim = mesh.topology.dim
    pyvista.start_xvfb()
    p = pyvista.Plotter()
    topology, cell_types, geometry = plot.create_vtk_mesh(mesh, tdim-1)
    marker = np.zeros(len(cell_types))
    marker[meshtags.indices] = meshtags.values
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid.cell_data["Marker"] = marker
    grid.set_active_scalars("Marker")
    p.add_mesh(grid, show_edges=True, cmap='jet', show_scalar_bar=False)
    p.show()


def plot_function(function):
    if not plot_in_line: return
    V = function.function_space
    p = pyvista.Plotter()
    topology, cell_types, geometry = plot.create_vtk_mesh(V)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    actor_0 = p.add_mesh(grid, style="wireframe", color="k")
    grid["u"] = function.x.array.reshape((geometry.shape[0], 3))
    warped = grid.warp_by_vector("u", factor=1.5)
    actor_1 = p.add_mesh(warped, show_edges=True)
    p.show()
