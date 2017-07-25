# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 17:29:36 2017

This script is intended to make nice plots of poincare spheres in the presence of PDL, in the spirit of Damask pp. 305.

@author: Noah
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch, ArrowStyle
from mpl_toolkits.mplot3d import proj3d

def add_arrows(axis, arrow_len):
    axis_arrow_length = arrow_len
    arrow = ArrowStyle("Fancy", 	head_length=0.4,head_width=0.2)
    arrow = 'simple'
    
    s1 = Arrow3D([-0.01, axis_arrow_length], [0, 0], [0, 0], mutation_scale=20,
                lw=1, arrowstyle=arrow, color="k")
    s2 = Arrow3D([0, 0], [-0.01, axis_arrow_length], [0, 0], mutation_scale=20,
                lw=1, arrowstyle=arrow, color="k")
    s3 = Arrow3D([0, 0], [0, 0], [-0.01, axis_arrow_length], mutation_scale=20,
                lw=1, arrowstyle=arrow, color="k")
    axis.add_artist(s1)
    axis.add_artist(s2)
    axis.add_artist(s3)


fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')
ax.set_axis_off()

u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:20j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x, y, z, color="k")
add_arrows(ax)
ax.set_ylim()


alpha = 0.345
direction = [-1, 0, 0]

dot_product = x*direction[0]*np.ones(np.shape(x)) + y*direction[1]*np.ones(np.shape(y)) +z*direction[2]*np.ones(np.shape(x)) 
trans = 1/(1+np.tanh(alpha)) * (1+np.tanh(alpha)*dot_product)

ax = fig.add_subplot(222, projection='3d')
ax.set_axis_off()
ax.plot_wireframe(trans*x, trans*y, trans*z, color="k")
add_arrows(ax)


class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args,                                   **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d,                renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)


