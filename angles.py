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

#top3left3 is the best tetrahedron design and top4left1 is the best principal polarization state design
if 0:
    design='top3left3'
    title='Design: Principal polarizations'
else:
    design='top4left1'
    title='Design: Regular Tetrahedron'
    
save_fig=1

if 'darwin' or 'linux' in sys.platform:
    directory = 'acquisition/data/small_metasurfaces/angles/' + design #top3left3'#top4left1' 
else:
    directory = 'acquisition\\data\\small_metasurfaces\\angles\\top4left1'

subdirs = ['order-2','order-1', 'order1', 'order2']
#angledirs = ['-40dgr','-38dgr','-36dgr','-34dgr','-32dgr','-30dgr','-28dgr','-26dgr','-24dgr','-22dgr','-20dgr','-18dgr','-16dgr','-14dgr','-12dgr','-10dgr','-8dgr','-6dgr','-4dgr','-2dgr','0dgr', '2dgr','4dgr','6dgr','8dgr','10dgr','12dgr','14dgr','16dgr','18dgr','20dgr','22dgr','24dgr','26dgr','28dgr','30dgr','32dgr','34dgr','36dgr','38dgr','40dgr']
angledirs = ['-36dgr','-32dgr','-28dgr','-24dgr','-20dgr','-16dgr','-12dgr','-8dgr','-4dgr',
             '0dgr','4dgr','8dgr','12dgr','16dgr','20dgr','24dgr','28dgr','32dgr','36dgr','40dgr']
os.chdir(directory)

data=[]
data_thorlabs=np.zeros((len(angledirs),len(subdirs),4))
indeks=0
for a in angledirs:
    index=0
    os.chdir(a)
#    print(os.path.dirname(os.path.realpath(__file__)))
    for folder in subdirs:
        os.chdir(folder)
        with open('polarimeter.txt') as f:
            polarimeter = np.array(list(csv.reader(f)), dtype='float')
            polarimeter = np.array([ 1 ] + list(np.mean(polarimeter,0)[0:3]))
            data_thorlabs[indeks][index][0:4]=polarimeter
        os.chdir('..') 
        index += 1
            

#    os.chdir('..')
    os.chdir('..')
    indeks+=1
    
    
data = np.array(data)    
data_thorlabs = np.array(data_thorlabs)


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
                antialiased=True, alpha=0.1, lw=0.)#, facecolors=cm)
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
   
t = np.linspace(0, 1, 21)  
farve = [np.array([1/3, 0, 0]),np.array([0, 1/3, 0]),np.array([0, 0, 1/3]),np.array([1/3, 1/3, 0])]
pp=[]
for i in range(4):   
    # Plotting Thorlabs polarimeter data
    S1 = data_thorlabs[:,i,1]
    S2 = data_thorlabs[:,i,2]
    S3 = data_thorlabs[:,i,3]
    
    #for i in range(0,4):
    #    for j in range(0,4):
    #        plt.plot(list(S1[n] for n in [i,j]),
    #                 list(S2[n] for n in [i,j]),
    #                 list(S3[n] for n in [i,j]), color='orange', lw=0.5, marker=' ')
             
    for j in range(len(S1)):
        pp.append(ax.plot([S1.item(j)/np.linalg.norm([S1.item(j),S2.item(j),S3.item(j)])], [S2.item(j)/np.linalg.norm([S1.item(j),S2.item(j),S3.item(j)])], [S3.item(j)/np.linalg.norm([S1.item(j),S2.item(j),S3.item(j)])],  c = farve[i]+2*j/(3*len(S1)), marker='o'))
#    ax.scatter(S1, S2, S3, marker='o')

#legend
#cbar = plt.colorbar()
#cbar.ax.set_yticklabels(['0','1','2','>3'])
#cbar.set_label('# of something', rotation=270)
red_proxy = plt.Rectangle((0, 0), 1, 1, fc=[2/3, 0, 0])
green_proxy = plt.Rectangle((0, 0), 1, 1, fc=[0, 2/3, 0])
blue_proxy = plt.Rectangle((0, 0), 1, 1, fc=[0, 0, 2/3])
yellow_proxy = plt.Rectangle((0, 0), 1, 1, fc=[2/3, 2/3, 0])
ax.legend([red_proxy,green_proxy,blue_proxy,yellow_proxy],['$m=-2$', '$m=-1$', '$m=+1$','$m=+2$'])
#plt.figlegend( (pp[0], pp[1], pp[2],pp[3]),
#    ('order-2', 'order-1', 'order1','order2'),
#    'upper right' )
ax.set_axis_off()
#plt.title(title)


if save_fig:
    file_name = 'StokesVectorsOnPSphere_'+ design + '.svg'
    os.chdir('../../../../../Graphics/angle')
    plt.savefig(file_name, format='svg')
#    os.chdir('..\\' + data_dir + '/' + linear_pol_extension)


plt.show()
       
helpline=[-1,1,-1,1]
for i in range(4):    
    S1 = data_thorlabs[:,i,1]
    S2 = data_thorlabs[:,i,2]
    S3 = data_thorlabs[:,i,3]
    # Turn off the axis planes
    #vinkler = np.arange(-40,41,2)
    vinkler = np.array([int(angle[:-3]) for angle in angledirs])
    #vinkler = np.insert(vinkler,0, [-40,-36,-32,-28,-24])
    #vinkler = np.insert(vinkler,len(vinkler), [24,28,32,36,40])
    plt.figure()
    #plt.plot(vinkler, S1,'-o')
    #plt.plot(vinkler, S2,'-o')
    #plt.plot(vinkler, S3,'-o')
    plot1=plt.scatter(vinkler, S1, color='blue')
    plot2=plt.scatter(vinkler, S2, color='orange')
    plot3=plt.scatter(vinkler, S3, color = 'green')
    plt.axhline(helpline[i], color='black', alpha=0.25)
    plt.figlegend( (plot1, plot2, plot3),
           ('$s_1$', '$s_2$', '$s_3$'),
           'upper right' )
    plt.title(subdirs[i])
    if save_fig:
        file_name = 'StokesVectorsVsAngle2D_'+ design + '_' + subdirs[i] + '.svg'
        plt.savefig(file_name, format='svg')
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
if ('left1' in directory) or ('left2' in directory):
    designed=np.array([[1,2*np.sqrt(2)/3,0,1/3],
                       [1,-np.sqrt(2)/3,np.sqrt(2/3),1/3],
                       [1,-np.sqrt(2)/3,-np.sqrt(2/3),1/3],
                       [1,0,0,-1]])
    twopsi_d=np.array([np.pi, 2*np.pi/3, -2*np.pi/3, 0])
    twochi_d=np.array([np.pi/6,np.pi/6, np.pi/6, -np.pi/2]) 
if ('left3' in directory) or ('left4' in directory):
    designed=np.array([[1,0,-1,0.01],
                       [1,0,0,1],
                       [1,0,0,-1],
                       [1,0,1,0.01]])
    twopsi_d=np.array([-np.pi/2, 0, 0, np.pi/2])
    twochi_d=np.array([0,np.pi/2, -np.pi/2, 0])

colrgb = np.array([0.5, 0.5, 0.5])    
n=0
p=[]
for q in range(len(angledirs)):

    data_thorlabs1=data_thorlabs[q,:,:]
    # iterate over measurement schemes
    for meas in [data_thorlabs1,designed]:
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
            if meas.all==data_thorlabs1.all:
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
            if meas.all==data_thorlabs1.all:
                linestyle='-'
                col=farve[1]+2*q/(3*len(angledirs))
                p.append(axarr[i].plot(xx,yy,color=col,ls=linestyle,alpha=0.8, label='Thorlabs measurement', linewidth=lw)[0])
            elif meas.all==designed.all:
                linestyle = '--'
                col = 'blue'
                p.append(axarr[i].plot(xx,yy, ls=linestyle, color=col, alpha=0.8, label='Designed states', linewidth=lw)[0])
                
            axarr[i].set_xlim([-1.1,1.1])
            axarr[i].set_ylim([-1.1,1.1])
            
            axarr[i].axes.get_xaxis().set_visible(False)
            axarr[i].axes.get_yaxis().set_visible(False)
    
            # find ang
    
            if abs(el) > 0.05 and (meas.all==data_thorlabs1.all or meas.all==designed.all):
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


plt.figlegend([p[0],p[7]],
           ['Thorlabs measurement','Designed states'],
           loc='lower center')


#axarr[1][0].axis('off')
#axarr[1][1].axis('off')
#axarr[1][2].axis('off')
#axarr[1][3].axis('off')
#plt.colorbar()
if save_fig:
    file_name = 'PolarizationEllipse_'+ design + '.svg'
    plt.savefig(file_name, format='svg')

plt.show()
