# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:42:59 2017

This is a script to analyze the angle dependence of the analyzer matrix

@contributors: Noah, Ruoping, Michael
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
from scipy import stats as stats
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.tri as mtri

plt.close('all')

save_fig = 0

linear_pol_extension = 'polarizer_only'  # folder for linear pol data
qwp_R = 'qwp_R'  # folder for qwp at first configuration
qwp_L = 'qwp_L'  # folder for qwp at second configuration

partial_pol = 'partial_pol8'  # folder location of partial pol data

plotStokes = 0 #do you want to plot Stokes parameters in cartesion (1) or polar (0) coordinates?

if 'linux' or 'darwin' in sys.platform:
    data_dir = 'acquisition/data/incident_angles_calibration'
else:
    data_dir = 'acquisition\\data\\incident_angles_calibration\\20deg'

os.chdir(data_dir)
#angles of (big) metasurface away from normal incidence
angledirs = ['0deg','5deg','10deg','15deg','20deg']

#%% Collect some error analysis functions
def covS(i, j, D, I, Dcov, Icov):
    ''' This function returns the element (i,j) of the covariance matrix of the result of Ainv*I
        D: the inverse instrument matrix
        I: the measured intensity vector
        Dcov: the covariance tensor for D
        Icov: the covariance matrix for I.
    '''
    assert len(I) == 4  # the measured intensity must have only four elements
    assert D.shape == (4, 4)  # the inverse instrument matrix is accordingly 4x4
    # the instrument matrix has is 4x4 so the covariance will be 4x4x4x4
    # because it is the covariance between two matrix elements, it requires four
    # inputs to fully specify each of which can range from 1-4     
    assert Dcov.shape == (4, 4, 4, 4)
    s = 0.0  # running total for covariance
    
    # add the first term from Eq. 22 in Ramos and Collados 
    for a in range(4):
        for b in range(4):            
            s += I[a] * I[b] * Dcov[i][a][j][b]
    # add the second term
    for k in range(4):
        for p in range(4):
            s += D[i][k] * D[j][p] * Icov[k][p]
    return s


def qwp_err(pd_arr):
    '''function returning list of qwp measurement averages for error analysis

    pd_arr: 1xn array with photodiode voltages, spanning the entire 360 degrees
    '''
    # must be an even # of measurements in the end and we divide total number by two                                     so it must be an even multiple of 4
    assert (len(pd_arr)/4) % 1 < 1e-12    
    avg_90deg = np.zeros(len(pd_arr))  # a vector to store the averages    
    offset = int(len(pd_arr)/4)  # number of steps to +90 deg measurement
    
    for i in range(len(pd_arr)):
        # % is so that measurements wrap around appropriately
        avg_90deg[i] = 0.5*(pd_arr[i]+pd_arr[(i+offset) % len(pd_arr)])
    
    return avg_90deg
    
def sorted_nicely( l ):
    """ From Mark Byers. Sorts the given iterable in the way that is expected.
 
    Required arguments:
    l -- The iterable to be sorted.
 
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key = alphanum_key)

def determine_dop(pol_state):
    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]

def errIncS(A_actual, A_perceived, Sinc):
    '''Calculate the measured Stokes vector (Sout) and error between Sinc and Sout

    Sout: measured Stokes vector (simulated)
    normSout: normalized Sout
    diffS: distance between Sinc and Sout vectors
    DiffDOP, diffS1, diffS2, diffS3: distance on each of 4 axes
    '''
    Imeas = np.dot(A_actual,Sinc)
    Sout = np.dot(A_perceived,Imeas)
    
    GreatCircledist = 1 #Choose between Eucledian 4D dist or Great Circle dist
    
    #euclidean distance plus great circle distance
    normSout = [Sout[1]/np.linalg.norm([Sout[1],Sout[2],Sout[3]]),Sout[2]/np.linalg.norm([Sout[1],Sout[2],Sout[3]]),Sout[3]/np.linalg.norm([Sout[1],Sout[2],Sout[3]])]
            #diffS = np.sqrt(np.square(Sinc[1]-normSout[0])+np.square(Sinc[2]-normSout[1])+np.square(Sinc[3]-normSout[2])) #Eucledean 3D method
    if GreatCircledist:
        diffS = np.arctan2(np.linalg.norm(np.cross(normSout,Sinc[1:])), np.dot(normSout,Sinc[1:]))
    else:
        #euclidean distance 4dim
        diffS = np.sqrt(np.square(Sinc[0]-Sout[0])+np.square(Sinc[1]-Sout[1])+np.square(Sinc[2]-Sout[2])+np.square(Sinc[3]-Sout[3]))
    
    
    SincDOP = determine_dop(Sinc)
    SoutDOP = determine_dop(Sout)
    
    
    diffDOP = np.abs(SincDOP-SoutDOP)
    diffS0 = np.abs(Sinc[0]-Sout[0])/Sinc[0]
    diffS1 = np.abs(Sinc[1]-normSout[0])
    diffS2 = np.abs(Sinc[2]-normSout[1])
    diffS3 = np.abs(Sinc[3]-normSout[2])
    
    azimutOut = 0.5 * np.arctan2(Sout[2],Sout[1])
    elliptOut = 0.5 * np.arctan2(Sout[3],np.sqrt(Sout[1]**2+Sout[2]**2))
    elliptOut2nd = Sout[3]/(Sout[0]+np.sqrt(Sout[1]**2+Sout[2]**2))
    
    azimutInc = 0.5 * np.arctan2(Sinc[2],Sinc[1])
    elliptInc = 0.5 * np.arctan2(Sinc[3],np.sqrt(Sinc[1]**2+Sinc[2]**2))
    elliptInc2nd = Sout[3]/(Sinc[0]+np.sqrt(Sinc[1]**2+Sinc[2]**2))
    
    diffAz = np.abs(azimutOut - azimutInc) #np.abs(np.mod(azimutOut - azimutInc + np.pi, 2*np.pi) - np.pi)
    diffEl = np.abs(elliptOut - elliptInc)
    diffA = azimutOut - azimutInc
    diffE = elliptOut - elliptInc
    
    #desperate actions to remove errors when crossing zero radians
    if diffAz>2.5 and Sout[1]<0 and np.abs(Sout[2])<0.1:
        diffAz=np.abs(diffAz-np.pi)
        diffA=diffAz
        

    # just difference
    #diffS=Sout-Sinc
    return normSout, diffS, diffDOP, diffS0, diffS1, diffS2, diffS3, diffAz, diffEl, diffA, diffE, Sout, GreatCircledist
#%% Short version of calibration.py
q=0
A=np.zeros((len(angledirs),4,4))
Ainv=np.zeros((len(angledirs),4,4))
#iterate over angles
for folder in angledirs:
    os.chdir(folder)

    angles = []  # angles off the linear polarizer
    inc_powers = []  # incident powers during measurement
    pd1_voltage = []  # photodiode 1 voltages
    pd2_voltage = []
    pd3_voltage = []
    pd4_voltage = []
    fit_functions = []
    fit_parameters = []
    variances = []
    
    os.chdir(linear_pol_extension)  # go get the linear pol data
    
    # iterate through the files in the directory which have been 'nicely' sorted
    for file in sorted_nicely(os.listdir()):
        if fnmatch.fnmatch(file, '*.txt'):
            params = file.split('_')
            try:  # if file name matches format do stuff with it
                power = float(params[-1][:-4])  # get rid of '.txt'
                angle = float(params[0][:-3])  # get rid of 'deg' in string
                inc_powers.append(power)
                angles.append(angle)
                # make file into array
                my_data = np.genfromtxt(file, delimiter=',')            
                pd1_voltage.append(np.mean(my_data[:, 0]))
                pd2_voltage.append(np.mean(my_data[:, 1]))
                pd3_voltage.append(np.mean(my_data[:, 2]))
                pd4_voltage.append(np.mean(my_data[:, 3]))
                
            except ValueError: # don't do anything with invalid file name
                pass
           
    # rearrange all data as sorted by pol_angle
    sorted_lists = sorted(zip(angles, inc_powers, pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage))
    # now recover each individual list
    angles, inc_powers, pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage =  [[x[i] for x in sorted_lists] for i in range(6)]
            
    # begin plotting of obtained data
    num_angles = len(angles) # number of angular test points
    angles = angles[:(num_angles)//2] # cut out part of array from 180 to 360 deg
    
    # split the photodiode powers into two separate lists
    
    pd1_voltage1 = np.array(pd1_voltage[:(num_angles)//2])
    pd2_voltage1 = np.array(pd2_voltage[:(num_angles)//2])
    pd3_voltage1 = np.array(pd3_voltage[:(num_angles)//2])
    pd4_voltage1 = np.array(pd4_voltage[:(num_angles)//2])
    pd1_voltage2 = np.array(pd1_voltage[(num_angles)//2:])
    pd2_voltage2 = np.array(pd2_voltage[(num_angles)//2:])
    pd3_voltage2 = np.array(pd3_voltage[(num_angles)//2:])
    pd4_voltage2 = np.array(pd4_voltage[(num_angles)//2:])
        
    # slice the incident power list in half as well
    inc_powers1 = np.array(inc_powers[:(num_angles)//2])
    inc_powers2 = np.array(inc_powers[(num_angles)//2:])
    
    # normalize each by the power incident during measurement
    pd1_voltage1 = pd1_voltage1/inc_powers1
    pd2_voltage1 = pd2_voltage1/inc_powers1
    pd3_voltage1 = pd3_voltage1/inc_powers1
    pd4_voltage1 = pd4_voltage1/inc_powers1
    pd1_voltage2 = pd1_voltage2/inc_powers2
    pd2_voltage2 = pd2_voltage2/inc_powers2
    pd3_voltage2 = pd3_voltage2/inc_powers2
    pd4_voltage2 = pd4_voltage2/inc_powers2
            
    # now average the two and propagate the error through the averaging
    pd1_voltage = (pd1_voltage1+pd1_voltage2)/2
    pd2_voltage = (pd2_voltage1+pd2_voltage2)/2
    pd3_voltage = (pd3_voltage1+pd3_voltage2)/2
    pd4_voltage = (pd4_voltage1+pd4_voltage2)/2
    
    # construct a larger matrix to hold the voltages
    pd_voltages = np.vstack((pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage))
    # convert these lists to proper numpy arrays
    angles = np.array(angles)
    inc_powers = np.array(inc_powers)
    
    
    # Now plot linear polarization cal data
    
    # define a function to plot the resulting graph
    def linear_cal_fig(xdata, ydata, min_angle, max_angle):
        
        # define a fitting function
        def fit_function(theta, a, b, c):
            return a + b*np.cos(2*theta*np.pi/180) + c*np.sin(2*theta*np.pi/180)
        
        fit_errs=[]
        fit_functions=[]
        # for photodiodes 1-4, fit data and store fitting parameters
        for i in range(0,4):
            x = xdata
            y = ydata[i, :]
            popt, variance = curve_fit(fit_function, x, y)  # fit a curve
    
            #standard error, assuming that the fit parameters are uncorrelated between each other
            standard_err = np.sqrt(np.diag(variance))
            fit_errs.append(standard_err)
            
            # store the fits as anonymous functions
            fit_functions.append(lambda theta: fit_function(theta, *popt))
            fit_parameters.append(popt)
            variances.append(variance)
    
        fit_errs=np.array(fit_errs)
        
        # now clip data at min and max angles
        mask = np.where((xdata>=min_angle) & (xdata<=max_angle))
        xdata = xdata[mask]
        ydata = ydata[:, mask]
            
#        for i in range(0, 4):
#            y = np.transpose(ydata[i, :])
        return fit_errs   
    
    save_figure = 0
    
    # now calc fit error
    min_angle = 0
    max_angle = 180
    fit_errs = linear_cal_fig(angles, pd_voltages, min_angle, max_angle)   
    
    # move onto the qwpR part of the calibration
    
    # get the data for both qwp sets of measurements
    for i in range(2):
        qwp_power_incR = []
        qwp_anglesR = []
        pol_anglesR = []
        pd1_voltageQR = []  # photodiode 1 voltages
        pd2_voltageQR = []
        pd3_voltageQR = []
        pd4_voltageQR = []
        
        os.chdir('..')  # move up a level
        if i == 1:  # determine which directory to cd into
            file = qwp_R
        else:
            file = qwp_L
        
        os.chdir(file)  # go get the qwp+pol data
        for file in sorted_nicely(os.listdir()):
            if fnmatch.fnmatch(file, '*.txt'):
                params = file.split('_')
                try:  # if file name matches format do stuff with it
                    power = float(params[-1][:-4])  # get rid of '.txt'
                    pol_angle = float(params[0][1:-3])  # get rid of 'deg' in string
                    qwp_angle = float(params[1][1:-3])
                    qwp_power_incR.append(power)
                    pol_anglesR.append(pol_angle)
                    qwp_anglesR.append(qwp_angle)
                    
                    # make file into array
                    data = np.genfromtxt(file, delimiter=',')
                    pd1_voltageQR.append(np.mean(data[:, 0]))
                    pd2_voltageQR.append(np.mean(data[:, 1]))
                    pd3_voltageQR.append(np.mean(data[:, 2]))
                    pd4_voltageQR.append(np.mean(data[:, 3]))
                except ValueError:  # don't do anything with invalid file name
                    print('file skipped ', file)
                    pass
        
        pol_anglesR = np.array(pol_anglesR)
        qwp_anglesR = np.array(qwp_anglesR)
        qwp_power_incR = np.array(qwp_power_incR)
        
        numAngles = len(pol_anglesR)
    
        # convert arrays to numpy format
        pd1_voltageQR = np.array(pd1_voltageQR)
        pd2_voltageQR = np.array(pd2_voltageQR)
        pd3_voltageQR = np.array(pd3_voltageQR)
        pd4_voltageQR = np.array(pd4_voltageQR)
        qwp_power_incR = np.array(qwp_power_incR)
        
        # normalize the photodiode voltages
        pd1_voltageQR = pd1_voltageQR/qwp_power_incR
        pd2_voltageQR = pd2_voltageQR/qwp_power_incR
        pd3_voltageQR = pd3_voltageQR/qwp_power_incR
        pd4_voltageQR = pd4_voltageQR/qwp_power_incR
 
        # average all the values in the list    
        if i == 0:
            # compute the mean value measured on each photodiode
            pd1R = np.mean(pd1_voltageQR)
            pd2R = np.mean(pd2_voltageQR)
            pd3R = np.mean(pd3_voltageQR)
            pd4R = np.mean(pd4_voltageQR)
            
        elif i == 1:
            pd1L = np.mean(pd1_voltageQR)
            pd2L = np.mean(pd2_voltageQR)
            pd3L = np.mean(pd3_voltageQR)
            pd4L = np.mean(pd4_voltageQR)
    
    # Construct the instrument matrix
    A3 = np.array([(pd1R-pd1L)/2, (pd2R-pd2L)/2, (pd3R-pd3L)/2, (pd4R-pd4L)/2])
    
    A3 = np.matrix.transpose(np.array([A3]))
    
    A02 = np.array(fit_parameters)  # data from linear calibration
    
    #Instrument matrix for each of q angles
    A[q][:][:] = np.hstack((A02, A3))  # This is the instrument matrix!
    
 #   print('Instrument matrix A: ')
 #   print(A)
 #   print('')
    
#    A_cond=np.linalg.cond(A[:][:][q], p = 2)  # condition number with the standard 2-norm
#    print('Condition number of A: ')
#    print(A_cond)
#    print('')

    #Inverse instrument matrix for each of q angles    
    Ainv[q][:][:] = np.linalg.inv(A[q][:][:])  # this is the inverse of the instrument matrix
#    print('Inverted matrix Ainv: ')
#    print(Ainv[:][:][q])
#    print('')

    
#    pickle.dump( Ainv_cov, open( "..\Ainv_cov.p", "wb" ) )
    os.chdir('..')
    os.chdir('..')
    q+=1

#plt.close('all')

#%% heatmap 

B=np.zeros((16,len(angledirs)-1))    
for indeks in range(1,len(angledirs)):
 #   print(indeks)
    taeller=0
    for i1 in range(4):
        for i2 in range(4):
            B[taeller][indeks-1]=abs((A[indeks][i1][i2]-A[0][i1][i2]))#/A[0][i1][i2])*100
            taeller+=1
plt.figure()
#z_min, z_max = 0,100# -np.abs(B).max(), np.abs(B).max()
plt.imshow(B)#,cmap='RdBu'extent=[5,20,17,1]
plt.xticks([0,1,2,3],('5','10','15','20'))
a=np.arange(15)
a=a[::2]
plt.yticks(a,('1','3','5','7','9','11','13','15'))
#plt.axis([0, 4, 0, 15])
#row_labels = list('WXYZ')
#ax.set_xticklabels(row_labels, minor=False)
plt.colorbar()
plt.xlabel('deg')
plt.ylabel('elements in A')
#plt.title('heatmap of deviation of A vs deg')

if save_fig:
    file_name = 'PolarimeterAngleHeatmapNotitle.svg'
    os.chdir('../../../Graphics/angle')
    plt.savefig(file_name, format='svg')

print('Analyzer matrix deviation: ')
print(B)

#%% Numerical calc of error on S

##OPTION 1: Symmetric sampling of sphere

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

##OPTION2: even sampling of sphere
#golden_angle = np.pi * (3 - np.sqrt(5)) 
#theta = golden_angle * np.arange(n) 
#z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
#radius = np.sqrt(1 - z * z)
# 
#x = np.zeros((n))
#t = np.zeros((n))
#x = radius * np.cos(theta)
#y = radius * np.sin(theta)

##OPTION3: spiral sampling of sphere
theta = np.linspace(0, 60*np.pi, n) 
if plotStokes:
    z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
else:
    z = np.linspace(0.95 - 1.0 / n, 1.0 / n - 0.95, n)
radius = np.sqrt(1 - z * z)
 
x = np.zeros((n))
y = np.zeros((n))
#t = np.zeros((n))
x = radius * np.cos(theta)
y = radius * np.sin(theta)

#delete values that arctan2 puts in wrong quadrant
#index=np.where(np.abs(y) < 0.05)
#x=np.delete(x, index)
#y=np.delete(y, index)
#z=np.delete(z, index)
#n=len(x)

#Define incoming Stokes vectors
Sinc=np.zeros((4,len(x)))
Sinc[0,:]=np.ones(len(x))
Sinc[1,:]=np.array(x)
Sinc[2,:]=np.array(y)
Sinc[3,:]=np.array(z)

#Calculate error on measured stokes vectors
S=[]
err=[]
dDOP=[]
dS0=[]
dS1=[]
dS2=[]
dS3=[]
dAz=[]
dEl=[]
dA=[]
dE=[]
Sout=[]
for j in range(1,len(angledirs)):
    for i in range(len(x)):
        SV, fejl, difdop, difS0, difS1, difS2, difS3, difAz, difEl, difA, difE, Sud, GreatCircledist  = errIncS(A[j][:][:], Ainv[0][:][:], Sinc[:,i])
        S.append(SV)
        err.append(fejl)
        dDOP.append(difdop)
        dS0.append(difS0)
        dS1.append(difS1)
        dS2.append(difS2)
        dS3.append(difS3)
        dAz.append(difAz)
        dEl.append(difEl)
        dA.append(difA)
        dE.append(difE)
        Sout.append(Sud)
S=np.array(S) #normalized measured Stokes
err=np.array(err) #error between measured and inc Stokes
S=np.transpose(S)
dS0=np.array(dS0) 
dS1=np.array(dS1) #error on S1
dS2=np.array(dS2) #error on S2
dS3=np.array(dS3) #error on S3
dAz=np.array(dAz) #abs error on Azimuthal angle
dEl=np.array(dEl) #abs error on ellipticity
dA=np.array(dA) #error on Azimuthal angle
dE=np.array(dE) #error on ellipticity
dDOP=np.array(dDOP) #error on SOP
Sout=np.array(Sout) #full measured Stokes vector



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

#plot incoming Stokes vectors without colormap
pSinc=ax.scatter3D(Sinc[1,:],Sinc[2,:],Sinc[3,:],color='red', marker='o',alpha=0.7)
ax.set_axis_off()
#lineInc = plt.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor='red',markersize=8)
plt.legend([pSinc],['Incoming polarization'])
if save_fig:
    file_name = 'PolarimeterAngleIncStokes.svg'
    plt.savefig(file_name, format='svg')
plt.show()

colbar=np.zeros((len(angledirs),len(x)))

def plot_color_sphere(ax, S, err, title, vinkel):
    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')

    plot_sphere(ax)
    maxred=np.max(err[(n*w):n*w+n]) #maximum error of meas with this angle
    mingreen=np.min(err[(n*w):n*w+n]) #minimum erro
    colbar[w,:]=err[(n*w):n*w+n]/maxred #normalize error to one
    
    #calculate root mean square error to compare
    if title=="Azimuthal error":
        rmse_dAz=np.sqrt(np.mean(err[(n*w):n*w+n]**2))
        rmse_dAz=np.around(rmse_dAz,decimals=4)
        #(mu, sigma) = stats.norm.fit(dA[(n*w):n*w+n])
        mu=np.mean(dA[(n*w):n*w+n])
        sigma=np.sqrt(np.mean((dA[(n*w):n*w+n]-mu)**2))
        sigma=np.around(sigma,decimals=4)
        title=title + ", $\sigma$=" + np.array2string(sigma) #np.array2string(rmse_dAz)
    #calculate root mean square error to compare    
    if title=="Ellipticity error":
        rmse_dEl=np.sqrt(np.mean(err[(n*w):n*w+n]**2))
        rmse_dEl=np.around(rmse_dEl,decimals=4)
        #(mu, sigma) = stats.norm.fit(dE[(n*w):n*w+n])
        mu=np.mean(dE[(n*w):n*w+n])
        sigma=np.sqrt(np.mean((dE[(n*w):n*w+n]-mu)**2))
        sigma=np.around(sigma,decimals=4)
        title=title + ", $\sigma$=" +  np.array2string(sigma)#np.array2string(rmse_dEl)
    
    farve=[colbar[w,:],1-colbar[w,:],np.zeros((len(colbar[w,:])))] #vary color from red (max error) to green (min error)
    # plot normalized measured stokes par (normSout) or incoming stokes par (Sinc)
    if plotMeasuredPol:
        ax.scatter3D([S[0,n*w:n*w+n]], [S[1,n*w:n*w+n]], [S[2,n*w:n*w+n]],c=np.transpose(farve), cmap='viridis',marker='o')#,alpha =0.5
    else:
        ax.scatter3D([Sinc[1,:]], [Sinc[2,:]], [Sinc[3,:]],c=np.transpose(farve), cmap='viridis',marker='o')#,alpha =0.5
    #plot perceived analyzer matrix
    for i in range(4):
        ax.scatter3D([A[0,i,1]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], [A[0,i,2]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], [A[0,i,3]/np.linalg.norm([A[0,i,1],A[0,i,2],A[0,i,3]])], color='cyan', marker='o')
    for i in range(0,4):
        for j in range(0,4):
            plt.plot(list(A[0,n,1]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]),
                     list(A[0,n,2]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]),
                     list(A[0,n,3]/np.linalg.norm([A[0,n,1],A[0,n,2],A[0,n,3]]) for n in [i,j]), color=[0.3,0.3,0.3], lw=1, marker=' ')
    
    #plot actual analyzer matrix 
    for i in range(4):
        ax.scatter3D([A[w+1,i,1]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], [A[w+1,i,2]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], [A[w+1,i,3]/np.linalg.norm([A[w+1,i,1],A[w+1,i,2],A[w+1,i,3]])], color='b', marker='o')
    for i in range(0,4):
        for j in range(0,4):
            plt.plot(list(A[w+1,n,1]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]),
                     list(A[w+1,n,2]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]),
                     list(A[w+1,n,3]/np.linalg.norm([A[w+1,n,1],A[w+1,n,2],A[w+1,n,3]]) for n in [i,j]), color='magenta', lw=1, marker=' ')

    ax.set_axis_off()
    #ax.set_title("max error (red): " + np.array_str(np.max(err[(800*w):800*w+800])))
    
    #make legend
    maxred=np.around(maxred, decimals=2)
    mingreen=np.around(mingreen, decimals=2)
    redtext="$max$ $error:$ " + np.array_str(maxred)
    greentext="$min$ $error:$ " + np.array_str(mingreen)
    line1 = plt.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor=[1,0,0],markersize=8, alpha=0.8)
    line2 = plt.Line2D(range(1), range(1), color="white", marker='o',markerfacecolor=[0,1,0],markersize=8, alpha=0.8)
    line3 = plt.Line2D(range(1), range(1), color="white", marker='o',markersize=8, markerfacecolor="blue")
    line4 = plt.Line2D(range(1), range(1), color="white", marker='o',markersize=8,markerfacecolor="cyan")
    yline = plt.Line2D(range(1), range(1), color="magenta")
    oline = plt.Line2D(range(1), range(1), color=[0.3,0.3,0.3])
    plt.legend((line1,line2,(line3, yline),(line4,oline)),(redtext,greentext, '$A_{actual}$', '$A_{perceived}$'),numpoints=1, loc="lower right")#A_actual', 'A_perceived'
    
    plt.title(title + ', angle = ' + vinkel)
    
    if save_fig:
        file_name = 'PolarimeterAngle' + vinkel + title + '.svg'
        plt.savefig(file_name, format='svg')
    plt.show()

plotMeasuredPol=1 #Choose whether you want to plot the incoming or (simulated) measured stokes vectors on the Pshere
#save_fig = 0
#os.chdir('../../../Graphics/angle/4D Eucledean dist')


#plot measured polarization and error shown as color
if GreatCircledist:
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, err,"Great Circle dist", angledirs[w+1]) #plot the Great Circle dist error
else:
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, err,"4D Eucledean dist", angledirs[w+1]) #or plot the Eucledean dist depending on choice in errIncS

if plotStokes:
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, dS1, "S1 error", angledirs[w+1]) # plot error on S1
        
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, dS2, "S2 error", angledirs[w+1]) # plot error on S2
        
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, dS3, "S3 error", angledirs[w+1]) # plot error on S3
else:
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, dAz, "Azimuthal error", angledirs[w+1]) # plot error on Azimuth
    for w in range(len(angledirs)-1):
        plot_color_sphere(plt.gca(), S, dEl, "Ellipticity error", angledirs[w+1]) # plot error on Azimuth

#for w in range(len(angledirs)-1):
#    plot_color_sphere(plt.gca(), S, dDOP, "DOP error", angledirs[w+1]) # plot error on DOP

#rmse_dAz=np.mean(dAz**2)
#rmse_dEl=np.mean(dEl**2)