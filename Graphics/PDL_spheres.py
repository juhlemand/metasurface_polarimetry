# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 17:29:36 2017

This script is intended to make nice plots of poincare spheres in the presence of PDL, in the spirit of Damask pp. 305.

@author: Noah
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch, ArrowStyle
from mpl_toolkits.mplot3d import proj3d

#%% Define a function for drawing the coordinate axes

arrow_len = 1.75

def add_arrows(axis, arrow_len):

    axis_arrow_length = arrow_len
    arrow = ArrowStyle("Fancy", 
              head_length=0.4, head_width=0.2)
    arrow = 'simple'
    s1 = Arrow3D([-0.01, axis_arrow_length], [0, 0], [0,
        0], mutation_scale=20, lw=1, arrowstyle=arrow,
        color="k")
    s2 = Arrow3D([0, 0], [-0.01, axis_arrow_length], [0,
        0], mutation_scale=20, lw=1, arrowstyle=arrow,
        color="k")
    s3 = Arrow3D([0, 0], [0, 0], [-0.01,
        axis_arrow_length], mutation_scale=20, lw=1,
        arrowstyle=arrow, color="k")
    
    axis.add_artist(s1)
    axis.add_artist(s2)
    axis.add_artist(s3)

#%% Set up the desired wire-frame plot. Return a handle to the axis in case we want to plot anything else there.

def wireframe_plot(x, y, z, subplot, figure):
    ax = figure.add_subplot(subplot, projection='3d')
    ax.set_axis_off()
    ax.plot_wireframe(x, y, z, color="k")
    add_arrows(ax, arrow_len)
    ax.set_xlim(-arrow_len/2, arrow_len/2)
    ax.set_ylim(-arrow_len/2, arrow_len/2)
    ax.set_zlim(-arrow_len/2, arrow_len/2)
    return ax

#%% Generate data for a normal Poincare Sphere

fig = plt.figure()  # open a figure to hold all future plots

# the initial Stokes vector
initial_state = np.array([0, -1, 0])
initial_state = initial_state/np.linalg.norm(initial_state)
arrow = ArrowStyle("Fancy", head_length=1, head_width=0.8)

init_state_vec = Arrow3D([0, initial_state[0]], [0, initial_state[1]], [0, initial_state[2]], lw=3, arrowstyle=arrow, color="r", mutation_scale=10)

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)

ax1 = wireframe_plot(x, y, z, 231, fig)  # plot the data
ax1.add_artist(init_state_vec)
#%% Now consider spherical plot of intensities under diattenuation
tx = 0.75
ty = 0.5

# construct the Jones matrix of the diattenuator
diattenuator = np.array([[tx, 0], [0, ty]])

# compute the attenuation factor
alpha = -np.log(ty/tx)  

direction = np.array([1, 0, 0])  # define the direction of diattentuation
gamma = np.tanh(alpha)

unit_matrix = np.ones(np.shape(x))

# compute a dot product between each vector and the direction, store as a matrix of dot products
dot_product = (x*direction[0] + y*direction[1]+        
               z*direction[2])

trans = 1/(1+gamma) * (1+gamma*dot_product)
initial_trans = 1/(1+gamma) * (1+gamma*np.dot(initial_state,  direction)) * initial_state

# now plot the result
ax2 = wireframe_plot(trans*x, trans*y, trans*z, 232, fig)
init_trans = Arrow3D([0, initial_trans[0]], [0, initial_trans[1]], [0, initial_trans[2]], lw=3, arrowstyle=arrow, color="r", mutation_scale=10)
ax2.add_artist(init_trans)

#%% Now compute the change in polarization state plot

# coefficients are long, so compute them separately
coeff_1 = np.sqrt(1-gamma**2)/(1+gamma*dot_product)
coeff_2 = gamma * (1+gamma**(-2)*(1-np.sqrt(1 - 
     gamma**2))*gamma*dot_product)/(1+gamma*dot_product)

new_x = coeff_1 * x + coeff_2 * direction[0]
new_y = coeff_1 * y + coeff_2 * direction[1]
new_z = coeff_1 * z + coeff_2 * direction[2]

# plot the result
ax3 = wireframe_plot(new_x, new_y, new_z, 233, fig)

dot_product = np.dot(initial_state, direction)
coeff_1 = np.sqrt(1-gamma**2)/(1+gamma*dot_product)
coeff_2 = gamma * (1+gamma**(-2)*(1-np.sqrt(1 - 
     gamma**2))*gamma*dot_product)/(1+gamma*dot_product) 

initial_trans2 = coeff_1*initial_state + coeff_2*direction
init_trans2 = Arrow3D([0, initial_trans2[0]], [0, initial_trans2[1]], [0, initial_trans2[2]], lw=3, arrowstyle=arrow, color="r", mutation_scale=10)
ax3.add_artist(init_trans2)

#%% Show a plot which shows polarization change, as well as the change in intensity

fin_x = trans*new_x
fin_y = trans*new_y
fin_z = trans*new_z

ax4 = wireframe_plot(fin_x, fin_y, fin_z, 234, fig)

initial_trans3 = 1/(1+gamma) * (1+gamma*np.dot(initial_state,  direction)) * initial_trans2
init_trans3 = Arrow3D([0, initial_trans3[0]], [0, initial_trans3[1]], [0, initial_trans3[2]], lw=3, arrowstyle=arrow, color="r", mutation_scale=10)
ax4.add_artist(init_trans3)
 

#%% Now simulate the effect of a retarder by rotating the resultant sphere

retardance = np.pi/2
ret_axis = [1, 0, 0] # retards about the x axis
ret_axis = ret_axis/np.linalg.norm(ret_axis)

shape_array = np.shape(fin_x)
rot_x = np.zeros(shape_array)
rot_y = np.zeros(shape_array)
rot_z = np.zeros(shape_array)

def rotate(axis, angle, vector):
    vect = np.array([curr_x, curr_y, curr_z])
    first_term = np.dot(axis, vector)*axis
    second_term = np.sin(angle)*np.cross(axis, vector)
    third_term = -np.cos(angle)*np.cross(axis, np.cross(axis, vector))
    return (first_term + second_term + third_term)


for i in range(shape_array[0]):
    for j in range(shape_array[1]):
      curr_x = fin_x[i, j]
      curr_y = fin_y[i, j]
      curr_z = fin_z[i, j]
      curr_vect = np.array([curr_x, curr_y, curr_z])
      
      final_vect = rotate(ret_axis, retardance, curr_vect)
      rot_x[i,j] = final_vect[0]
      rot_y[i,j] = final_vect[1]
      rot_z[i,j] = final_vect[2]

ax5 = wireframe_plot(rot_x, rot_y, rot_z, 235, fig)

initial_trans4 = rotate(ret_axis, retardance, initial_trans3)
init_trans4 = Arrow3D([0, initial_trans4[0]], [0, initial_trans4[1]], [0, initial_trans4[2]], lw=3, arrowstyle=arrow, color="r", mutation_scale=10)
ax5.add_artist(init_trans4)

#%% A class found on StackExchange for drawing 3D arrows in matplotlib

class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0),
                                 *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d,
                     zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)