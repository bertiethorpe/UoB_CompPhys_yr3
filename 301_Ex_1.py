#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 14:30:04 2022

@author: bertiethorpe
"""

#======# Computational Physics 301 21/22 - Exercise 1 - Linear Algebra #======#

# Program to find the tension in three cables suspending a skycam, such as seen
# in sports stadiums. This is done through solving 3 linear equations with the 
# scipy linalg.solve family of functions. 

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as linalg # Used over equivalent scipy.linalg for array broadcasting functionality

G = 9.81 # ms^-2
M = 50 # kg. mass of camera
Z = -7 # m. constant vertical distance of hanging camera
DELTA_Y = 0.5 # difference element

weight = np.array([[0] ,[0] ,[-M*G]])

# Set up array of cable anchor points on 2D plane. 
anchor = np.array([[-45*np.sqrt(3), -45], [0, 90], [45*np.sqrt(3), -45]])

# Initialise bounding coordinates from max/min values of the triangle
_, y_min = min(anchor, key = lambda m: m[1])
_, y_max = max(anchor, key = lambda m: m[1])
x_min, _  = min(anchor, key = lambda m: m[0])
x_max, _  = max(anchor, key = lambda m: m[0])

y_vals = np.arange(y_min, y_max+DELTA_Y, DELTA_Y)
x_vals = np.linspace(x_min, x_max, np.size(y_vals))

# Create meshgrid of coords, then unravel to 1D array for ease of computation
y_cam, x_cam = np.meshgrid(y_vals, x_vals, indexing='ij')
coords = np.column_stack((x_cam.ravel(), y_cam.ravel()))

def point_in_triangle():
    """
    Function to determine whether coordinate is within bounds of plane triangle,
    calculated using barycentric coordinate system. Returns 2D boolean array
    the same dimensions as the bounding coords.
    """
    # Initialise barycentric coords relative to the first anchor point
    vecB = anchor[1] - anchor[0]
    vecC = anchor[2] - anchor[0]
    vecP = coords - anchor[0]
    
    dotBB = np.dot(vecB, vecB)
    dotCC = np.dot(vecC, vecC)
    dotCB = np.dot(vecC, vecB)
    dotCP = np.dot(vecP, vecC)
    dotBP = np.dot(vecP, vecB)
    
    inv_denom = 1 / (dotCC * dotBB - dotCB * dotCB)
    
    # Calculate barycentric weight constants
    alpha = (dotBB * dotCP - dotCB * dotBP) * inv_denom
    beta = (dotCC * dotBP - dotCB * dotCP) * inv_denom
    total = alpha + beta
    
    # Compare weights to generate bool list
    array = (alpha >= 0) & (beta >= 0) & (total < 1)
    
    grid = np.reshape(array, (len(y_vals),len(x_vals)), 'F')
    
    return grid

boolean_grid = point_in_triangle()
print(boolean_grid)

# Mask camera meshgrid with the boolean array
x_masked = x_cam * boolean_grid 
y_masked = y_cam * boolean_grid

def cable_tension():
    """
    
    """
    # Angles of each cable
    theta_0 = np.arctan((y_masked - anchor[0,1]) / (x_masked - anchor[0,0]))
    theta_1 = np.arctan((x_masked - anchor[1,0]) / (anchor[1,1] - y_masked)) + np.radians(30)
    theta_2 = np.arctan((y_masked - anchor[2,1]) / (anchor[2,0] - x_masked))
    
    # Vertical angles 
    phi_0 = np.arctan(Z/np.sqrt((x_masked - anchor[0,0])**2 + (y_masked - anchor[0,1])**2))
    phi_1 = np.arctan(Z/np.sqrt((x_masked - anchor[1,0])**2 + (anchor[1,1] - y_masked)**2))
    phi_2 = np.arctan(Z/np.sqrt((anchor[2,0] - x_masked)**2 + (y_masked - anchor[2,1])**2))
    
    # Array of matrix form of system of equations 
    matrix = np.array([[-np.cos(theta_0), np.sin(theta_1 - np.radians(30)), np.cos(theta_2)],
                       [-np.sin(theta_0), np.cos(theta_1 - np.radians(30)), -np.sin(theta_2)],
                       [np.sin(phi_0), np.sin(phi_1), np.sin(phi_2)]
                       ])
    
    #matrix = np.nan_to_num(matrix) # filter out NaN matrix values
    matrix = np.transpose(matrix) # transpose so linalg.solve can iterate
    
    weight_tile = np.resize(weight,(271,271,3,1))
    T = linalg.solve(matrix, weight_tile)
    
    return T

tension = cable_tension()

print(tension)

fig, axs = plt.subplots(1, 3, figsize=(15, 3), sharey=True)
axs[0].pcolormesh(x_cam, y_cam, (tension[:,:,0,0] / np.rot90(boolean_grid, k=3)))
axs[1].pcolormesh(x_cam, y_cam, (tension[:,:,1,0] / np.rot90(boolean_grid, k=3)))
axs[2].pcolormesh(x_cam, y_cam, (tension[:,:,2,0] / np.rot90(boolean_grid, k=3)))




