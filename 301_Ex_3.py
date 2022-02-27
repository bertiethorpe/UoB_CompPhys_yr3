#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 18:23:52 2022

@author: bertiethorpe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import scipy.sparse as sparse
from scipy.sparse.linalg import spsolve

RHO = 7900 # density, kg m^-3
LAMBDA = 59 # thermal conductivity, W m^-1 K^-1
SHC = 450 # specific heat capacity, J kg^-1 K^-1

# define thermal diffusivity
therm_diff = LAMBDA / (RHO * SHC)

T_ROOM = 293.15
T_HOT = 1273.15
T_COLD = 273.15

LENGTH = 0.5 # of rod, m
N = 100 # number of discrete space points

# define space and time steps
delta_x = LENGTH / N
delta_t = 10

alpha = (therm_diff * delta_t) / (delta_x ** 2)

def dirichlet_interior_nodes():
    """ 
    Function returns matrix for unknown interior nodes of the 1D heat equation
    Dirichlet boundary problem.
    """
    # define the interior diagonals for the sparse matrix
    diagonals = [[np.full((N-2),-alpha)],
                  [np.full((N-2),(1 + 2*alpha))],
                  [np.full((N-2),-alpha)]]
    
    matrix = sparse.diags(diagonals, [-1,0,1], shape=(N-2,N-2))
    return matrix

def neumann_interior_nodes():
    """ 
    Function returns matrix for unknown interior nodes of the 1D heat equation
    with one end a Dirichlet boundary and the other a Neumann boundary.
    """
    # define the interior diagonals for the sparse matrix
    diagonals = [[np.full((N-2),-alpha)],
                  [np.full((N-2),(1 + 2*alpha))],
                  [np.full((N-2),-alpha)]]
    
    matrix = sparse.diags(diagonals, [-1,0,1], shape=(N-2,N-2))
    
    # change element to reflect neumann boundary
    # must use compressed row sparse matrix to ammend values
    csrmatrix =  matrix.tocsr()
    csrmatrix[-1,-1] = 1 + alpha
    matrix = csrmatrix.todia()
    
    return matrix

def dirichlet_boundary(T0,TN):
    """
    Function returns vector containing the Dirichlet boundary conditions
    for both ends of 1D rod.
    """
    vector = np.zeros((N-2))
    vector[0] = -alpha*T0
    vector[-1] = -alpha*TN
    return vector

def neumann_boundary(T0):
    """
    Function returns vector containing the Neumann boundary conditions for 
    one end of 1D rod. The other end resorts to a Dirichlet boundary."""
    vector = np.zeros((N-2))
    vector[0] = -alpha*T0
    return vector

dirbound = dirichlet_boundary(T_HOT,T_COLD)
neubound = neumann_boundary(T_HOT)

def heat_diffusion(matrix, boundary):
    
    # initialise temp values, simulation time, etc. for time iteration
    T = T_ROOM * np.ones(N-2)
    t_span = 3000
    rod = np.linspace(delta_x/2, LENGTH-delta_x/2, N-2)
    
    # setup array to store temps at each time step
    results = []
    
    fig, ax = plt.subplots()
    
    for n in np.arange(t_span/delta_t +1):
        
        # start by clearing the axes
        ax.cla()
        
        ax.set_xlim(0, LENGTH)
        ax.set_ylim(T_COLD, T_HOT)
        ax.set_xlabel('Rod Distance (m)')
        ax.set_ylabel('Temperature (K)')
        
        sol = spsolve(matrix,T-boundary)
        results.append(sol)
    
        # reinitialise temps with updated array
        T = sol
        
        # multi-colour line setup
        points = np.array([rod, T]).T.reshape(-1,1,2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(T_COLD,T_HOT)
        lc = LineCollection(segments, cmap='plasma',norm=norm)
        lc.set_array(T)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        plt.show()
        plt.pause(0.0001)
  

heat_diffusion(dirichlet_interior_nodes(),dirbound)
heat_diffusion(neumann_interior_nodes(),neubound)

#=================== DISCUSSION ========================#

# Scipy sparse diagonal matrix routines were used for efficient creation of 
# tridiagonal matrices. 

# Scipy.sparse.linalg.spsolve takes a sparse matrix to solve the  matrix equation 
# u_n = M * u_n+1 + b, obtained through a Taylor expansion solving to find the 
# temperature at each  finite difference position element along the rod. This 
# function was preferable to the general linalg.solve as the matrices are 
# specifically input and can be ensured to be sparse. If there were any special 
# cases, the more general solution would be preferable as spsolve executes very 
# slowly without sparse matrices.

# A simple time-stepping loop was employed to evolve the diffusion of heat along 
# the rod over time. 

# an OOP approach was taken given the two similar boundary cases.
# The first was a Dirichlet boundary set of functions, setting the ends of the 
# rod at constant temp. The second was a half Dirichlet, half von Neumann 
# boundaries; one of the ends was kept at zero flux out of the rod.

# The results from the Dirichlet diffusion showed to approach constant gradient 
# temp distribution from one end to the other. The results from the Neumann 
# diffusion showed a slow increase in temp across the whole rod; left running 
# long enough, the rod would approach a uniform temperature (that of the furnace
# end). This was expected given the nature of Neumann flux.

