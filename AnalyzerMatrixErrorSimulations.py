# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:42:59 2017

This is a script to analyze polarimetry calibration data.

@contributors: Noah, Ruoping
"""
#%reset -f

import os, re, pickle
import fnmatch
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

plt.close('all')

#linear_pol_extension = 'polarizer_only'  # folder for linear pol data
#qwp_R = 'qwp_R'  # folder for qwp at first configuration
#qwp_L = 'qwp_L'  # folder for qwp at second configuration
#
#partial_pol = 'partial_pol8'  # folder location of partial pol data
#
#power_meter_error = 0.001 #Error in power meter reading from ambient light, unit in mW
#
#
#if 'linux' or 'darwin' in sys.platform:
#    data_dir = 'acquisition/data/incident_angles_calibration'
#else:
#    data_dir = 'acquisition\\data\\incident_angles_calibration\\20deg'
#
#os.chdir(data_dir)
#angledirs = ['0deg','5deg','10deg','15deg','20deg']

#%% Collect some error analysis functions



def determine_dop(pol_state):
    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]

def errIncS(A_actual, A_perceived, Sinc):
    '''forklaring

    variabler
    '''
    Imeas = np.dot(A_actual,Sinc)
    Sout = np.dot(A_perceived,Imeas)
    #Sout=np.array(Sout)
    
    #euclidean distance 3dim plus great circle distance dim
    normSout = [Sout[1]/np.linalg.norm([Sout[1],Sout[2],Sout[3]]),Sout[2]/np.linalg.norm([Sout[1],Sout[2],Sout[3]]),Sout[3]/np.linalg.norm([Sout[1],Sout[2],Sout[3]])]
    #diffS = np.sqrt(np.square(Sinc[1]-normSout[0])+np.square(Sinc[2]-normSout[1])+np.square(Sinc[3]-normSout[2]))
    GSdiffS = np.arctan2(np.linalg.norm(np.cross(normSout,Sinc[1:])), np.dot(normSout,Sinc[1:]))
    
    SincDOP = determine_dop(Sinc)
    SoutDOP = determine_dop(Sout)
    
    diffDOP = np.abs(SincDOP-SoutDOP)
    diffS1 = np.abs(Sinc[1]-Sout[1])
    diffS2 = np.abs(Sinc[2]-Sout[2])
    diffS3 = np.abs(Sinc[3]-Sout[3])
    
    maxdiffS = GSdiffS+diffDOP
    
    #euclidean distance 4dim
    diffS = np.sqrt(np.square(Sinc[0]-Sout[0])+np.square(Sinc[1]-Sout[1])+np.square(Sinc[2]-Sout[2])+np.square(Sinc[3]-Sout[3]))

    # just difference
    #diffS=Sout-Sinc
    return normSout, diffS, diffDOP, diffS1, diffS2, diffS3, Sout, maxdiffS

#%% Define Analyzer matrices
m=2
A=np.zeros((m,4,4))
Ainv=np.zeros((m,4,4))
A[0][:][:]=np.array([[1,0,0,-1],
           [1,-np.sqrt(2)/3,-np.sqrt(2/3),1/3],
           [1,-np.sqrt(2)/3,np.sqrt(2/3),1/3],
           [1,2*np.sqrt(2)/3,0,1/3]])
A[1][:][:]=A[0][:][:]
A[1,0,0]=A[0,0,0]+0.1

##rotate around S1
#v=10*np.pi/180
#rotx=np.array([[1,0,0,0],
#           [0,1,0,0],
#           [0,0,np.cos(v),-np.sin(v)],
#           [0,0,np.sin(v),np.cos(v)]])
#A[1][:][:]=np.transpose(np.dot(rotx,np.transpose(A[0][:][:])))#np.dot(A[0][:][:],rotx)

#A[1][:][:]=A[0][:][:]
#A[1,:,0]=A[1,:,0]+0.1


for i in range(len(A)):
    Ainv[i][:][:]=np.linalg.inv(A[i][:][:])



#%% heatmap 
plt.close('all')
B=np.zeros((16,len(A)-1))    
for indeks in range(1,len(A)):
 #   print(indeks)
    taeller=0
    for i1 in range(4):
        for i2 in range(4):
            B[taeller][indeks-1]=abs((A[indeks][i1][i2]-A[0][i1][i2]))#/A[0][i1][i2])*100
            taeller+=1
plt.figure()
#z_min, z_max = 0,100# -np.abs(B).max(), np.abs(B).max()
plt.imshow(B)#,cmap='RdBu'
#plt.axis([0, 4, 0, 15])
#row_labels = list('WXYZ')
#ax.set_xticklabels(row_labels, minor=False)
plt.colorbar()
plt.xlabel('deg')
plt.ylabel('elements in A')
plt.title('heatmap of deviation (in %) of A vs deg')
#print('Analyzer matrix deviation: ')
#print(B)

#%% Numerical calc of error on S

##Symmetric sampling of sphere

#noofpphi=20
#noofpchi=40
#phiangle = np.linspace(0, 0.5*np.pi, noofpphi)
#chiangle = np.linspace(0, np.pi, noofpchi)
#
#n=noofpchi*noofpphi
#
##phiangle, chiangle = np.meshgrid(phiangle, chiangle)
#
## The Cartesian coordinates of the unit sphere
#x=[]
#y=[]
#z=[]
#for a in range(len(phiangle)):
#    for b in range(len(chiangle)):
#        x.append(np.cos(2*phiangle[a]) * np.cos(2*chiangle[b]))
#        y.append(np.sin(2*phiangle[a]) * np.cos(2*chiangle[b]))
#        z.append(np.sin(2*chiangle[b]))
        
n = 800

##even sampling of sphere
#golden_angle = np.pi * (3 - np.sqrt(5)) 
#theta = golden_angle * np.arange(n) 
#z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
#radius = np.sqrt(1 - z * z)
# 
#x = np.zeros((n))
#t = np.zeros((n))
#x = radius * np.cos(theta)
#y = radius * np.sin(theta)

#spiral sampling of sphere
theta = np.linspace(0, 60*np.pi, n) 
z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
radius = np.sqrt(1 - z * z)
 
x = np.zeros((n))
t = np.zeros((n))
x = radius * np.cos(theta)
y = radius * np.sin(theta)



Sinc=np.zeros((4,len(x)))
Sinc[0,:]=np.ones(len(x))
Sinc[1,:]=np.array(x)
Sinc[2,:]=np.array(y)
Sinc[3,:]=np.array(z)

#,,np.array(y),np.array(z)]
#Sinc=np.array(Sinc)
#Sinc=np.matrix.transpose(Sinc)
#Sforskel=np.zeros((len(phiangle),len(angledirs)))
S=[]
err=[]
dDOP=[]
dS1=[]
dS2=[]
dS3=[]
Sout=[]
maxfejl=[]
for j in range(1,len(A)):
    for i in range(len(x)):
        SV, fejl, difdop, difS1, difS2, difS3, Sud, maksfejl  = errIncS(A[j][:][:], Ainv[0][:][:], Sinc[:,i])
        S.append(SV)
        err.append(fejl)
        dDOP.append(difdop)
        dS1.append(difS1)
        dS2.append(difS2)
        dS3.append(difS3)
        Sout.append(Sud)
        maxfejl.append(maksfejl)
        #Sforskel[i][j]=errIncS(A[:][:][j], Ainv[:][:][0], Sinc[:,i])
S=np.array(S)
err=np.array(err)
S=np.transpose(S)        
dS1=np.array(dS1)
dS2=np.array(dS2)
dS3=np.array(dS3)
dDOP=np.array(dDOP)
Sout=np.array(Sout)
maxfejl=np.array(maxfejl)



#%%  plot
# plotting on Poincare sphere

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

# Set the aspect ratio to 1 so our sphere looks spherical
fig = plt.figure(figsize=plt.figaspect(1.))
ax = fig.add_subplot(111, projection='3d')

def plot_sphere(ax,arrows='xyz',equatorial=True):
    phi = np.linspace(0, np.pi, 200)
    theta = np.linspace(0, 2*np.pi, 200)

    #equatorial circle
    xe=np.sin(theta)
    ye=np.cos(theta)

    phi, theta = np.meshgrid(phi, theta)

    # The Cartesian coordinates of the unit sphere
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)

    ax.plot_surface(x, y, z,  rstride=10, cstride=10, color='#EBE3E8',
                antialiased=True, alpha=0.1, lw=0.5)#, facecolors=cm)
    if 'y' in arrows:
        ax.add_artist(Arrow3D([0, 0], [-0.03, 1.5], 
                        [0,0], mutation_scale=15, 
                        lw=0.25, arrowstyle="-|>", color="black"))
        ax.text(0,1.5,0, '$S_2$', fontweight='bold')        
    if 'x' in arrows:
        ax.add_artist(Arrow3D([0.0, 1.5], [0,0], 
                        [0,0], mutation_scale=15, 
                        lw=0.25, arrowstyle="-|>", color="black"))
        ax.text(1.6,0,0, '$S_1$', fontweight='bold')        
    if 'z' in arrows:        
        ax.add_artist(Arrow3D([0, 0], [0,0], 
                        [-0.03,1.5], mutation_scale=15, 
                        lw=0.25, arrowstyle="-|>", color="black"))
        ax.text(0,0,1.5, '$S_3$',fontweight='bold')
    if equatorial:
        ax.plot(xe,ye,0,'--', dashes=(10, 10), lw=0.25, color='red', alpha=1)
    
plot_sphere(ax)

#for j in range(len(x)):
pSinc=ax.scatter3D(Sinc[1,:],Sinc[2,:],Sinc[3,:],color='red', marker='o',alpha=0.7)
ax.set_axis_off()
#lineInc = plt.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor='red',markersize=8)
plt.legend([pSinc],['Incoming polarization'])
plt.show()

colbar=np.zeros((len(A),len(x)))

def plot_color_sphere(ax, S, err, title):
    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')

    plot_sphere(ax)
    maxred=np.max(err[(n*w):n*w+n])
    mingreen=np.min(err[(n*w):n*w+n])
    colbar[w,:]=err[(n*w):n*w+n]/maxred
    
    farve=[colbar[w,:],1-colbar[w,:],np.zeros((len(colbar[w,:])))]
    ax.scatter3D([S[0,n*w:n*w+n]], [S[1,n*w:n*w+n]], [S[2,n*w:n*w+n]],c=np.transpose(farve), cmap='viridis',marker='o')#,alpha =0.5
        
    for i in range(4):
        ax.scatter3D([A[0,i,1]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], [A[0,i,2]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], [A[0,i,3]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], color='cyan', marker='o')
        
    for i in range(0,4):
        for j in range(0,4):
            plt.plot(list(A[0,n,1]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]),
                     list(A[0,n,2]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]),
                     list(A[0,n,3]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]), color=[0.3,0.3,0.3], lw=1, marker=' ')
    
    
    for i in range(4):
        ax.scatter3D([A[w+1,i,1]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], [A[w+1,i,2]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], [A[w+1,i,3]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], color='b', marker='o')

    for i in range(0,4):
        for j in range(0,4):
            plt.plot(list(A[w+1,n,1]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]),
                     list(A[w+1,n,2]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]),
                     list(A[w+1,n,3]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]), color='magenta', lw=1, marker=' ')

    ax.set_axis_off()
    #ax.set_title("max error (red): " + np.array_str(np.max(err[(800*w):800*w+800])))
    maxred=np.around(maxred, decimals=3)
    mingreen=np.around(mingreen, decimals=3)
    
    redtext="max error: " + np.array_str(maxred)
    greentext="min error: " + np.array_str(mingreen)
    
#    red_proxy = plt.Rectangle((0, 0), 1, 1, fc=[1, 0, 0])
#    green_proxy = plt.Rectangle((0, 0), 1, 1, fc=[0, 1, 0])
#    blue_proxy = plt.Rectangle((0, 0), 1, 1, fc=[0, 0, 1])
#    black_proxy = plt.Rectangle((0, 0), 1, 1, fc=[0, 0, 0])
#    ax.legend([red_proxy,green_proxy,blue_proxy,black_proxy],[redtext, greentext, 'A_actual (yellow tetrahedron)','A_perceived (orange tetrahedron)'],loc='lower center')

    line1 = plt.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=[1,0,0],markersize=8)
    line2 = plt.Line2D(range(1), range(1), color="white", marker='o',markerfacecolor=[0,1,0],markersize=8)
    line3 = plt.Line2D(range(1), range(1), color="white", marker='o',markersize=8, markerfacecolor="blue")
    line4 = plt.Line2D(range(1), range(1), color="white", marker='o',markersize=8,markerfacecolor="cyan")
    yline = plt.Line2D(range(1), range(1), color="magenta")
    oline = plt.Line2D(range(1), range(1), color=[0.3,0.3,0.3])
    plt.legend((line1,line2,(line3, yline),(line4,oline)),(redtext,greentext, 'A_actual', 'A_perceived'),numpoints=1, loc=8)
    plt.title(title)
    plt.show()


for w in range(len(A)-1):
    plot_color_sphere(plt.gca(), S, err,"Great circle error, change one analyzer vector")
    
for w in range(len(A)-1):
    plot_color_sphere(plt.gca(), S, dS1, "S1 error, change one analyzer vector")
    
for w in range(len(A)-1):
    plot_color_sphere(plt.gca(), S, dS2, "S2 error, change one analyzer vector")
    
for w in range(len(A)-1):
    plot_color_sphere(plt.gca(), S, dS3, "S3 error, change one analyzer vector")

for w in range(len(A)-1):
    plot_color_sphere(plt.gca(), S, dDOP, "DOP error, change one analyzer vector")

