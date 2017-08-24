import matplotlib.pyplot as plt
import numpy as np
import sys, csv, os


from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.colors import LightSource
import random

if 'linux' or 'darwin' in sys.platform:
    directory = 'acquisition/data/big_metasurface'
    directory = 'acquisition/data/small_metasurfaces/top2_left3'
else:
    directory = 'acquisition\\data\\small_metasurfaces\\top5_left4'

subdirs = ['order_-2', 'order_-1', 'order_1', 'order_2']
os.chdir(directory)

####################################################################
#calculating efficiency
try:
    with open('powers.txt') as f:
        raw = f.readlines()

    efficiency = 0.
        
    for i in range(len(raw)):
        raw[i] = raw[i].split(' ')
        if 'mW' in raw[i][-1]:
            raw[i][1] = float(raw[i][1])*1000
        if 'inc' in raw[i][0]:
            inc_power = float(raw[i][1])
        if raw[i][0]!='0:' and raw[i][0]!='inc:':
            efficiency += float(raw[i][1])
            
    efficiency = efficiency/inc_power
    print('4-orders efficiency:',efficiency)
except:
    pass

####################################################################
#%% calculating efficiency

data=[]
data_thorlabs=[]

for folder in subdirs:
    os.chdir(folder)
    # if time sequential measurement was done, fetch the data
    try:
        with open('pol_only') as f:
            pol_only = np.array(list(csv.reader(f)), dtype='float')
        with open('qwp_L') as f:
            qwp_L = np.array(list(csv.reader(f)), dtype='float')
            qwp_L = np.mean(qwp_L,0)
        with open('qwp_R') as f:
            qwp_R = np.array(list(csv.reader(f)), dtype='float')
            qwp_R = np.mean(qwp_R,0)
    except:
        pass
    
    # get data from ThorLabs polarimeter
    with open('polarimeter.txt') as f:
        polarimeter = np.array(list(csv.reader(f)), dtype='float')
        polarimeter = np.array([ 1 ] + list(np.mean(polarimeter,0)[0:3]))
        data_thorlabs.append(polarimeter)

    try:
        data.append(np.hstack([pol_only.transpose()[1],qwp_R.transpose()[1],qwp_L.transpose()[1]]))
    except:
        pass
    os.chdir('..')
    
    
data = np.array(data)    
data_thorlabs = np.array(data_thorlabs)

A = [[1,1,0,0],
     [1,0,1,0],
     [1,-1,0,0],
     [1,0,-1,0],
     [1,0,0,1],
     [1,0,0,-1]]
A = np.array(A)
Ainv = np.linalg.pinv(A)

# apply the pseudo-inverse to the data measured in time sequential measurement
measured_stokes=np.zeros((4,4))
for i in range(len(data)):
   measured_stokes[i] = np.dot(Ainv, data[i])
   measured_stokes[i][0] = 1.
   measured_stokes[i][1:] = measured_stokes[i][1:] / np.linalg.norm(measured_stokes[i][1:])


#############################################################################
#%% plotting on Poincare sphere

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
                antialiased=True, alpha=0.5, lw=0.)#, facecolors=cm)
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

# Plotting Thorlabs polarimeter data
S1 = data_thorlabs.transpose()[1]
S2 = data_thorlabs.transpose()[2]
S3 = data_thorlabs.transpose()[3]
print(np.vstack([S1,S2,S3]).transpose())

# iterate over all pairs of Stokes coordinates
for i in range(0,4):
    for j in range(0,4):
        plt.plot(list(S1[n] for n in [i,j]),
                 list(S2[n] for n in [i,j]),
                 list(S3[n] for n in [i,j]), color='orange', lw=0.5, marker=' ')

# Plotting measured polarization data
S1 = measured_stokes.transpose()[1]
S2 = measured_stokes.transpose()[2]
S3 = measured_stokes.transpose()[3]

# iterate over all pairs of Stokes coordinates
for i in range(0,4):
    for j in range(0,4):
        plt.plot(list(S1[n] for n in [i,j]),
                 list(S2[n] for n in [i,j]),
                 list(S3[n] for n in [i,j]), color='blue', lw=0.5, marker=' ')

# Turn off the axis planes
ax.set_axis_off()
plt.show()

#########################################################
#%% Plotting polarization ellipses
fig_x = 15
fig_y = fig_x/3
my_dpi=192
fig, axarr = plt.subplots(1,4, figsize=(fig_x, fig_y), dpi=my_dpi)

if 'big_metasurface' in directory:
    designed=np.array([[1,0,0,-1],
                       [1,-np.sqrt(2)/3,-np.sqrt(2/3),1/3],
                       [1,-np.sqrt(2)/3,np.sqrt(2/3),1/3],
                       [1,2*np.sqrt(2)/3,0,1/3]])
    twopsi_d=np.array([0, -2*np.pi/3, 2*np.pi/3, np.pi])
    twochi_d=np.array([-np.pi/2, np.pi/6,np.pi/6, np.pi/6])
if 'left1' or 'left2' in directory:
    designed=np.array([[1,2*np.sqrt(2)/3,0,1/3],
                       [1,-np.sqrt(2)/3,np.sqrt(2/3),1/3],
                       [1,-np.sqrt(2)/3,-np.sqrt(2/3),1/3],
                       [1,0,0,-1]])
    twopsi_d=np.array([np.pi, 2*np.pi/3, -2*np.pi/3, 0])
    twochi_d=np.array([np.pi/6,np.pi/6, np.pi/6, -np.pi/2]) 
if 'left3' or 'left4' in directory:
    designed=np.array([[1,0,-1,0.01],
                       [1,0,0,1],
                       [1,0,0,-1],
                       [1,0,1,0.01]])
    twopsi_d=np.array([-np.pi/2, 0, 0, np.pi/2])
    twochi_d=np.array([0,np.pi/2, -np.pi/2, 0])
    
n=0
p=[]

# iterate over measurement schemes
for meas in [data_thorlabs,measured_stokes,designed]:
    S0 = meas.transpose()[0]
    S1 = meas.transpose()[1]
    S2 = meas.transpose()[2]
    S3 = meas.transpose()[3]

    #S0 = np.array([1,1])
    #S1 = np.array([0,0])
    #S2 = np.array([0,0])
    #S3 = np.array([1,-1])

    S1=S1/S0
    S2=S2/S0
    S3=-S3/S0
    err=0.
    # iterate over diffraction orders
    for i in range(4):
        twopsi = np.arctan2(S2[i], S1[i])
        twochi = np.arctan2(S3[i], np.sqrt(S1[i]**2+S2[i]**2))

        def diff(a,b):
            d1=abs(a+b)
            d2=abs(a-b)
            d3=abs(a+b+np.pi)
            d4=abs(a+b-np.pi)
            return min([d1,d2,d3,d4])
        
        # if we are in the Thorlabs data case
        if meas.all==data_thorlabs.all:
            err += diff(twopsi, twopsi_d[i])**2+diff(twochi,twochi_d[i])**2
        # find the ellipticity, the tangent of the ellipticity angle
        el=np.tan((twochi)/2)
        
        # set direction in which ellipse is traversed
        if np.sign(S3[i]) > 0:
            thetas = np.linspace(0, 2*np.pi, 10000)
        else:
            thetas = np.linspace(2*np.pi, 0, 10000)
            
        b = 1 # assume major axis is 1
        a = b/np.abs(el) # now compute the length of the minor axis
        norm = max(a,b) # normalize by whichever was larger
        a=a/norm
        b=b/norm
        # I think this results in being off by 90 degrees:        
        #r = a*b/np.sqrt((a*np.cos(thetas)**2+(b*np.sin(thetas))**2))

        #x= r*np.cos(thetas)
        #y= r*np.sin(thetas)
        
        x = a*np.cos(thetas)
        y = b*np.sin(thetas)        
        
        # rotate the ellipse
        xx=x*np.cos(0.5*twopsi)-y*np.sin(0.5*twopsi)
        yy=x*np.sin(0.5*twopsi)+y*np.cos(0.5*twopsi)
        
        axarr[i].set_aspect('equal')        
        
        lw = 1.5
        # plot the polarization ellipse
        if meas.all==data_thorlabs.all:
            linestyle='-'
            col='red'
            p.append(axarr[i].plot(xx,yy,color=col,ls=linestyle,alpha=0.8, label='Thorlabs measurement', linewidth=lw)[0])
        elif meas.all==measured_stokes.all:
            linestyle = '--'
            col = 'red'
            p.append(axarr[i].plot(xx,yy,color=col, ls=linestyle, alpha=0.8, label='Time-sequential measurement', linewidth=lw)[0])
        elif meas.all==designed.all:
            linestyle = '--'
            col = 'blue'
            p.append(axarr[i].plot(xx,yy, ls=linestyle, color=col, alpha=0.8, label='Designed states', linewidth=lw)[0])
            
        axarr[i].set_xlim([-1.1,1.1])
        axarr[i].set_ylim([-1.1,1.1])
        
        axarr[i].axes.get_xaxis().set_visible(False)
        axarr[i].axes.get_yaxis().set_visible(False)

        # find ang

        if abs(el) > 0.05 and (meas.all==data_thorlabs.all or meas.all==designed.all):
            points = [2000, 4500, 8000]
            u = (np.roll(xx, 1) - xx)[points]
            v = (np.roll(yy, 1) - yy)[points]
            norms = np.sqrt(u**2+v**2)
            u, v = u/norms, v/norms # equalize arrow lengths            
            
            axarr[i].quiver(xx[points], yy[points], u, v, color=col, headwidth=10, headlength=10, pivot='mid', headaxislength=10)            
            
#            axarr[0][i].arrow(xx[3*len(xx)//8-10], yy[3*len(yy)//8-10],
#                           np.sign(S3[i])*(xx[3*len(xx)//8+10]-xx[3*len(xx)//8]),
#                           np.sign(S3[i])*(yy[3*len(yy)//8+10]-yy[3*len(xx)//8]),
#                           head_width=0.1, head_length=0.2, linewidth=0., alpha=.8, color=col)
#            axarr[0][i].arrow(xx[7*(len(xx)//8-10)], yy[7*(len(yy)//8-10)],
#                           np.sign(S3[i])*(xx[7*(len(xx)//8)+10]-xx[7*(len(xx)//8)]),
#                           np.sign(S3[i])*(yy[7*(len(yy)//8)+10]-yy[7*(len(xx)//8)]),
#                           head_width=0.1, head_length=0.2, linewidth=0., alpha=.8, color=col)
    n+=1
    if err>0.00000001:
        print('Error with respect to designed', err)


plt.figlegend([p[0],p[4],p[8]],
           ['Thorlabs measurement','Time-sequential measurement','Designed states'],
           loc='lower center')


#axarr[1][0].axis('off')
#axarr[1][1].axis('off')
#axarr[1][2].axis('off')
#axarr[1][3].axis('off')
plt.show()



###########################################################################
#%% Prepare figure for export
left = os.getcwd() # directory to come back to
os.chdir('../../../..')
os.chdir('Graphics')
names = directory.split('\\')
new_dir = names[-1]
try:
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    os.chdir(new_dir)
    f_name = 'ellipses.pdf'
    plt.savefig(f_name, dpi=my_dpi)
except:
    pass
os.chdir(left)
