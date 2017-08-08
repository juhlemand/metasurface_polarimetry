import csv, os, pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from sys import platform
from scipy.stats import circmean
from scipy import stats as stats

polarimeter_file = 'polarimeter.txt' #polarimeter data file
DOP_CUTOFF=0.5

#instrument matrix from calibration
if 'linux' in platform:
    directory='acquisition/data/calibration4/comparison' #data location folder
    os.chdir(directory)
    Ainv=np.loadtxt('../Ainv.txt')
    Ainv_cov=pickle.load(open( "../Ainv_cov.p", "rb" ))
else:
    directory='acquisition\\data\\calibration4\\comparison' #data location folder
    os.chdir(directory)
    Ainv=np.loadtxt('..\\Ainv.txt')
    Ainv_cov=pickle.load(open( "..\\Ainv_cov.p", "rb" ))

N_measurements=len(os.listdir())-1 #number of measurements which have been taken

#https://www.ruhr-uni-bochum.de/ika/forschung/forschungsbereich_kolossa/Daten/Buchkapitel_Uncertainty.pdf

def covS(i,j, D, I, Dcov, Icov):
    ''' This function returns the element (i,j) of the covariance matrix of the result of Ainv*I
        D: the inverse instrument matrix
        I: the measured intensity vector
        Dcov: the covariance tensor for D
        Icov: the covariance matrix for I.
    '''
    assert len(I)==4
    assert D.shape==(4,4)
    assert Dcov.shape==(4,4,4,4)
    s=0.0
    for a in range(4):
        for b in range(4):            
            s += I[a]*I[b]*Dcov[i][a][j][b]
    for k in range(4):
        for l in range(4):
            s += D[i][k]*D[j][l]*Icov[k][l]
    return s

class polarimeter:
    def __init__(self, N_measurements):
        self.N_measurements = N_measurements
        self.data = np.zeros((N_measurements,4))
        self.stdev = np.zeros((N_measurements,4))
        self.poincare = np.zeros((N_measurements,3))
        self.poincare_std = np.zeros((N_measurements,3))
        
    def add(self, index, stokes_vector, stokes_stdev, poincare, poincare_std):
        assert stokes_vector.shape ==(4,)
        assert stokes_stdev.shape ==(4,)
        self.data[index] = stokes_vector
        self.stdev[index] = stokes_stdev
        self.poincare[index] = poincare
        self.poincare_std[index] = poincare_std

    def __str__(self):
        print(self.data)
        return str(self.stdev)
        
polarimeter_raw = []
with open(polarimeter_file, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter = ',')
    for row in reader:
        #print(', '.join(row))
        polarimeter_raw.append(row)
polarimeter_raw = np.array(polarimeter_raw)

#####################################################
#%% Parsing polarimeter file
polarimeter_data = polarimeter(N_measurements)
p = []
is_measurement = 0
index = 0

for row in polarimeter_raw:
    # measurements begin with this start flag in text file
    if row[:15] == ['#####START#####']:
        is_measurement = 1
    # measurements end with this end flag in the text file
    elif row[:13] == ['#####END#####']:
        is_measurement = 0
        
        # these are tasks to carry out once all the data in a given file has been extracted
        s_v=np.array([1./(0.01*np.mean(np.array(p),0)[-2]), np.mean(np.array(p),0)[0],np.mean(np.array(p),0)[1],np.mean(np.array(p),0)[2]])
        s_std=np.array([0.01*np.std(np.array(p),0)[-2]/(s_v[0])**2, np.std(np.array(p),0)[0],np.std(np.array(p),0)[1],np.std(np.array(p),0)[2]])
        poincare=np.array([0.01*np.mean(np.array(p),0)[-2],np.mean(np.array(p),0)[3],np.mean(np.array(p),0)[4]])
        poincare_std=np.array([0.01*np.std(np.array(p),0)[-2],np.std(np.array(p),0)[3],np.std(np.array(p),0)[4]])
        polarimeter_data.add(index, s_v, s_std, poincare, poincare_std)
        index += 1
        p=[]
        
    elif is_measurement:
        p.append(np.array(row,dtype=float))

###############################################################
#%% Parsing metasurface data
fnames = os.listdir()
fnames.remove(polarimeter_file)
fnames = sorted(fnames, key=lambda item: (int(item.partition('_')[0])
                               if item[0].isdigit() else float('inf'), item))

if len(fnames) != polarimeter_data.N_measurements:
    raise ValueError

metasurface_data, cov_m,= [], []
cov_m = []
for n in range(len(fnames)):
    temp = []
    with open(fnames[n], 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            temp.append([])
            for el in row:
                temp[-1].append(float(el))
    temp = np.array(temp)
    metasurface_data.append(np.average(temp,0))
    I = np.average(temp, 0)
    # assume covariance matrix is diagonal
    Icov = np.diag(np.std(temp, 0)**2)
    #covariance of stokes vector
    cov_stokes = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            cov_stokes[i][j]=covS(i, j, Ainv, I, Ainv_cov, Icov)
    #taking the diagonal as standard error of stokes vector
    cov_m.append(cov_stokes)
        
cov_m = np.array(cov_m)  # cov_m is a list of covariance matrices
metasurface_data = np.array(metasurface_data)

###############################################################
# Analyzing and plotting comparison

for i in range(len(metasurface_data)):
    metasurface_data[i] = np.dot(Ainv, metasurface_data[i].transpose())

m_dops = np.sqrt(metasurface_data.transpose()[1]**2+metasurface_data.transpose()[2]**2+metasurface_data.transpose()[3]**2)/metasurface_data.transpose()[0]
p_dops = np.sqrt(polarimeter_data.data.transpose()[1]**2+polarimeter_data.data.transpose()[2]**2+polarimeter_data.data.transpose()[3]**2)/polarimeter_data.data.transpose()[0]

#error in dop
# dop = sqrt(S1**2+S2**2+S3**2)/S0

S3 = polarimeter_data.data.transpose()[3]
S2 = polarimeter_data.data.transpose()[2]
S1 = polarimeter_data.data.transpose()[1]
S0 = polarimeter_data.data.transpose()[0]
dS3 = polarimeter_data.stdev.transpose()[3]
dS2 = polarimeter_data.stdev.transpose()[2]
dS1 = polarimeter_data.stdev.transpose()[1]
dS0 = polarimeter_data.stdev.transpose()[0]
p_dops_err = np.sqrt((dS0*np.sqrt(S1**2+S2**2+S3**2)/S0**2)**2
                   +(dS1*S1/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2
                   +(dS2*S2/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2
                   +(dS3*S3/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2)

#metasurface
S3 = metasurface_data.transpose()[3]
S2 = metasurface_data.transpose()[2]
S1 = metasurface_data.transpose()[1]
S0 = metasurface_data.transpose()[0]
dS3 = np.sqrt(cov_m[:,3,3])
dS2 = np.sqrt(cov_m[:,2,2])
dS1 = np.sqrt(cov_m[:,1,1])
dS0 = np.sqrt(cov_m[:,0,0])
cov_S0_S1 = cov_m[:,0,1]
cov_S0_S2 = cov_m[:,0,2]
cov_S0_S3 = cov_m[:,0,3]
cov_S2_S1 = cov_m[:,2,1]
cov_S3_S1 = cov_m[:,3,1]
cov_S3_S2 = cov_m[:,3,2]
#https://www.wolframalpha.com/input/?i=derivative+of+sqrt(x1**2%2Bx2**2%2Bx3**2)%2Fx0
m_dops_err=np.sqrt((dS0*np.sqrt(S1**2+S2**2+S3**2)/S0**2)**2
                   +(dS1*S1/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2
                   +(dS2*S2/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2
                   +(dS3*S3/(S0*np.sqrt(S1**2+S2**2+S3**2)))**2
                   +2*cov_S0_S1*(-S1/S0**3)
                   +2*cov_S0_S2*(-S2/S0**3)
                   +2*cov_S0_S3*(-S3/S0**3)
                   +2*cov_S2_S1*(S2*S1/(S0**2 * (S1**2+S2**2+S3**2)))
                   +2*cov_S3_S1*(S3*S1/(S0**2 * (S1**2+S2**2+S3**2)))
                   +2*cov_S3_S2*(S3*S2/(S0**2 * (S1**2+S2**2+S3**2)))
                   )


plt.errorbar(range(0,len(fnames)),m_dops, alpha=0.5, yerr=m_dops_err, label='metasurface',fmt='.')
plt.errorbar(range(0,len(fnames)),p_dops, alpha=0.5, yerr=p_dops_err, label='thorlabs',fmt='.')
plt.plot((0,len(fnames)), (1.0,1.0), color='gray', alpha=0.5)
plt.ylim([0,1.1])
plt.legend()
plt.show()

##############################################################
#plotting code

# histogram parameters
histogram_color = (0,1,1, 0.5)
fit_line_width = 3
legend_size = 12

# inset plotting parameters
xmin_dop = 0.5
xmax_dop = 0.6
ymin_dop = 0.5
ymax_dop = 0.6

xmin_azimuth = 1.5
xmax_azimuth = 2
ymin_azimuth = 1.5
ymax_azimuth = 2

xmax_altitude = 0
xmin_altitude = 0.25
ymin_altitude = 0
ymax_altitude = 0.25

edge_color = "0.0"
zoom = 4
location = 4 # lower right corner


#poincare sphere coordinates in radians
# azimuth psi
m_2psi=np.arctan2(S2, S1)
p_2psi=np.arctan2(polarimeter_data.data.transpose()[2], polarimeter_data.data.transpose()[1])
#p_2psi=np.pi*polarimeter_data.poincare.transpose()[1]/180
#stdev error in S2/S1:
#http://123.physics.ucdavis.edu/week_9_files/taylor_209-226.pdf
p_2psi_err=2*np.pi*polarimeter_data.poincare_std.transpose()[1]/180
m_2psi_err=np.sqrt((S1*dS2/(S1**2+S2**2))**2
                   +(S2*dS1/(S1**2+S2**2))**2
                   -2*S1*S2*cov_S2_S1/(S1**2+S2**2)**2)

# altitude chi
#p_2chi=np.arctan(polarimeter_data.data.transpose()[3]/np.sqrt(polarimeter_data.data.transpose()[1]**2+polarimeter_data.data.transpose()[2]**2))
p_2chi=2*np.pi*polarimeter_data.poincare.transpose()[2]/180
m_2chi=-np.arctan2(S3, np.sqrt(S1**2+S2**2))
#https://www.wolframalpha.com/input/?i=derivative+of+atan(x3%2Fsqrt(x1**2%2Bx2**2))
m_2chi_err=np.sqrt((dS1*S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS2*S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS3*np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))**2
                   +2*cov_S2_S1*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S1*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S2*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   )
p_2chi_err=2*np.pi*polarimeter_data.poincare_std.transpose()[2]/180

#removing from DOP_CUTOFF
for i in range(len(m_dops)):
    if m_dops[i]<DOP_CUTOFF and p_dops[i]<DOP_CUTOFF:
        for arr in [m_dops,p_dops,m_dops_err,p_dops_err,
                    m_2psi,p_2psi,m_2psi_err,p_2psi_err,
                    m_2chi,p_2chi,m_2chi_err,p_2chi_err]:
            arr[i]=9e999-9e999 #nan

data_arr=[m_dops,p_dops,m_dops_err,p_dops_err,
          m_2psi,p_2psi,m_2psi_err,p_2psi_err,
          m_2chi,p_2chi,m_2chi_err,p_2chi_err]
for i in range(len(data_arr)):
    data_arr[i]=data_arr[i][~np.isnan(data_arr[i])]
m_dops=data_arr[0]
p_dops=data_arr[1]
m_dops_err=data_arr[2]
p_dops_err=data_arr[3]
m_2psi=data_arr[4]
p_2psi=data_arr[5]
m_2psi_err=data_arr[6]
p_2psi_err=data_arr[7]
m_2chi=data_arr[8]
p_2chi=data_arr[9]
m_2chi_err=data_arr[10]
p_2chi_err=data_arr[11]

#fixing small line segment
for i in range(len(m_2psi)):
    if m_2psi[i] < -2 and p_2psi[i] > 2:
        m_2psi[i]=m_2psi[i]+2*np.pi

# shift points due to azimuth offset
diffs = m_2psi - p_2psi
p_2psi = p_2psi + circmean(diffs, high=np.pi, low=-np.pi)
p_2psi = p_2psi + 2*np.pi*(p_2psi < -np.pi)
# now correct positions of weird outlying points - using * as elementwise and
m_2psi = m_2psi + 2*np.pi*((m_2psi<0)*(p_2psi>0))
m_2psi = m_2psi + 2*np.pi*((m_2psi>0)*(p_2psi<0))
# recalculate diffs for histogram
diffs = np.mod(m_2psi - p_2psi + np.pi, 2*np.pi) - np.pi
# now remove points whose error bars are inordinately large
mask = np.logical_and((p_2psi_err < 3*np.mean(p_2psi_err)), (m_2psi_err <    3*np.mean(m_2psi_err)))
diffs = diffs[mask]
m_2psi = m_2psi[mask]
p_2psi = p_2psi[mask]
m_2psi_err = m_2psi_err[mask]
p_2psi_err = p_2psi_err[mask]
m_2chi = m_2chi[mask]
p_2chi = p_2chi[mask]
m_2chi_err = m_2chi_err[mask]
p_2chi_err = p_2chi_err[mask]
m_dops = m_dops[mask]
p_dops = p_dops[mask]
m_dops_err = m_dops_err[mask]
p_dops_err = p_dops_err[mask]

# having cleaned up the data, we plot the results for dop, chi, and psi
# plot the dop data
f, axarr  = plt.subplots(2,3)
axarr[0][0].scatter(p_dops, m_dops,alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].errorbar(p_dops, m_dops, xerr=p_dops_err, yerr=m_dops_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].plot([0,1],[0,1],alpha=0.75,color='black')
axarr[0][0].set_title('DOP')
axarr[0][0].set_xlabel('Polarimeter measurement')
axarr[0][0].set_ylabel('Metasurface measurement')
axarr[0][0].set_xlim([-0.1,1.1])
axarr[0][0].set_ylim([-0.1,1.1])

# create a zoomed inset of the data
axins_dop = zoomed_inset_axes(axarr[0][0], zoom, loc=location)
axins_dop.scatter(p_dops, m_dops, alpha=0.5, s=2.)
axins_dop.errorbar(p_dops, m_dops, xerr=p_dops_err, yerr=m_dops_err, alpha=0.5, fmt=' ')
axins_dop.plot([-np.pi,np.pi],[-np.pi, np.pi], color='black', alpha=0.75)
x1, x2, y1, y2 = xmin_dop, xmax_dop, ymin_dop, ymax_dop
axins_dop.set_xlim(x1, x2) # apply the x-limits
axins_dop.set_ylim(y1, y2) # apply the y-limits
mark_inset(axarr[0][0], axins_dop, loc1=3, loc2=1, fc="none", ec=edge_color)
axins_dop.set_xticklabels([])
axins_dop.set_yticklabels([])

# plot a histogram
diffs=m_dops-p_dops
n, bins, patches = axarr[1][0].hist(diffs, bins=np.linspace(min(diffs), max(diffs) + 0.005, np.sqrt(len(diffs))), facecolor = histogram_color)
axarr[1][0].axvline(0.0,color='black', alpha=0.25)
(mu, sigma) = stats.norm.fit(diffs)
sample = np.linspace(min(bins), max(bins), 200)
y = mlab.normpdf(sample, mu, sigma)
l = axarr[1][0].plot(sample, len(diffs)*(bins[1]-bins[0])*y, 'r--', linewidth=fit_line_width, label= '$\mu=%.3f$\n$\sigma=%.3f $'%(mu, sigma))
axarr[1][0].legend(prop={'size': legend_size})

# move on to azimuth

axarr[0][1].scatter(p_2psi, m_2psi, alpha=0.5, s=2.)
axarr[0][1].errorbar(p_2psi, m_2psi, xerr=p_2psi_err, yerr=m_2psi_err, alpha=0.5, fmt=' ')
axarr[0][1].set_title('Azimuth $2\chi$')
axarr[0][1].set_xlabel('Polarimeter measurement (radians)')
axarr[0][1].set_ylabel('Metasurface measurement (radians)')
axarr[0][1].set_xlim([-1.1*np.pi,1.1*np.pi])
axarr[0][1].set_ylim([-1.1*np.pi,1.1*np.pi])
axarr[0][1].plot([-np.pi,np.pi],[-np.pi, np.pi], color='black', alpha=0.75)

# create a zoomed inset of the data
axins_az = zoomed_inset_axes(axarr[0][1], zoom, loc=location)
axins_az.scatter(p_2psi, m_2psi, alpha=0.5, s=2.)
axins_az.errorbar(p_2psi, m_2psi, xerr=p_2psi_err, yerr=m_2psi_err, alpha=0.5, fmt=' ')
axins_az.plot([-np.pi,np.pi],[-np.pi, np.pi], color='black', alpha=0.75)
x1, x2, y1, y2 = xmin_azimuth, xmax_azimuth, ymin_azimuth, ymax_azimuth
axins_az.set_xlim(x1, x2) # apply the x-limits
axins_az.set_ylim(y1, y2) # apply the y-limits
mark_inset(axarr[0][1], axins_az, loc1=2, loc2=1, fc="none", ec=edge_color)
axins_az.set_xticklabels([])
axins_az.set_yticklabels([])

# histogram
n, bins, patches = axarr[1][1].hist(diffs, bins=np.linspace(min(diffs), max(diffs) + 0.005, np.sqrt(len(diffs))), facecolor=histogram_color)
axarr[1][1].axvline(0.0,color='black', alpha=0.25)
(mu, sigma) = stats.norm.fit(diffs)
sample = np.linspace(min(bins), max(bins), 200)
y = mlab.normpdf(sample, mu, sigma)
l = axarr[1][1].plot(sample, len(diffs)*(bins[1]-bins[0])*y, 'r--', linewidth=fit_line_width, label= '$\mu=%.3f$\n$\sigma=%.3f $'%(mu, sigma))
axarr[1][1].legend(prop={'size': legend_size})
axarr[1][1].axvline(0.0,color='black', alpha=0.25)

# move on to altitude

popt = np.polyfit(p_2chi, m_2chi, 1)
if popt[0]<0: # correct for RCP/LCP ambiguity
    m_2chi=-m_2chi
    
axarr[0][2].scatter(p_2chi, m_2chi, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].errorbar(p_2chi, m_2chi, xerr=p_2chi_err, yerr=m_2chi_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].set_title('Altitude $\phi$')
axarr[0][2].set_xlabel('Polarimeter measurement (radians)')
axarr[0][2].set_ylabel('Metasurface measurement (radians)')
axarr[0][2].set_xlim([-1.1*np.pi/2,1.1*np.pi/2])
axarr[0][2].set_ylim([-1.1*np.pi/2,1.1*np.pi/2])
axarr[0][2].plot([-np.pi/2,np.pi/2],[-np.pi/2, np.pi/2], color='black', alpha=0.75)

# create a zoomed inset of the data
axins_alt = zoomed_inset_axes(axarr[0][2], zoom, loc=location)
axins_alt.scatter(p_2chi, m_2chi, alpha=0.5, s=2.)
axins_alt.errorbar(p_2chi, m_2chi, xerr=p_2chi_err, yerr=m_2chi_err, alpha=0.5, fmt=' ')
axins_alt.plot([-np.pi/2,np.pi/2],[-np.pi/2, np.pi/2], color='black', alpha=0.75)
x1, x2, y1, y2 = xmin_altitude, xmax_altitude, ymin_altitude, ymax_altitude
axins_alt.set_xlim(x1, x2) # apply the x-limits
axins_alt.set_ylim(y1, y2) # apply the y-limits
mark_inset(axarr[0][2], axins_alt, loc1=4, loc2=2, fc="none", ec=edge_color)
axins_alt.set_xticklabels([])
axins_alt.set_yticklabels([])
axins_alt.tick_params(axis='y', direction='in')

# plot a histogram
diffs=m_2chi-p_2chi
n, bins, patches = axarr[1][2].hist(diffs, bins=np.linspace(min(diffs), max(diffs) + 0.005, np.sqrt(len(diffs))), facecolor = histogram_color)
axarr[1][2].axvline(0.0,color='black', alpha=0.25)
(mu, sigma) = stats.norm.fit(diffs)
sample = np.linspace(min(bins), max(bins), 200)
y = mlab.normpdf(sample, mu, sigma)
l = axarr[1][2].plot(sample, len(diffs)*(bins[1]-bins[0])*y, 'r--', linewidth=fit_line_width, label= '$\mu=%.3f$\n$\sigma=%.3f $'%(mu, sigma))
axarr[1][2].legend(prop={'size': legend_size})
axarr[1][2].axvline(0.0,color='black',alpha=0.25)

plt.show()

savefig = 0
f_name = 'comparison_graphs.svg'
if savefig:
    plt.savefig(f_name)


#%% Poincare sphere plots
#############################################################################
# plotting on Poincare sphere

from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib.colors import LightSource
import random

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

def plot_sphere(ax,arrows='xyz', equatorial=True):
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

    ax.plot_surface(x, y, z,  rstride=2, cstride=2, color='#EBE3E8',
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
# Plotting selected datapoints
npoints=50
dpoints=[]
for n in range(npoints):
    dpoints.append(int(random.random()*len(m_dops)))

dpoints=np.array(dpoints)

for n in range(npoints):
    S3=np.sin(m_2chi[dpoints[n]])
    S2=np.sin(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]])
    S1=np.cos(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]])
    ax.add_artist(Arrow3D([0, S1/np.linalg.norm([S1,S2,S3])],
                          [0, S2/np.linalg.norm([S1,S2,S3])],
                          [0, S3/np.linalg.norm([S1,S2,S3])], mutation_scale=5,
                          lw=0.5, arrowstyle="-|>", color="blue"))

    x=np.linspace(np.cos(m_2psi[dpoints[n]]-m_2psi_err[dpoints[n]])*np.cos(m_2chi[dpoints[n]]),
                  np.cos(m_2psi[dpoints[n]]+m_2psi_err[dpoints[n]])*np.cos(m_2chi[dpoints[n]]))
    y=np.linspace(np.sin(m_2psi[dpoints[n]]-m_2psi_err[dpoints[n]])*np.cos(m_2chi[dpoints[n]]),
                  np.sin(m_2psi[dpoints[n]]+m_2psi_err[dpoints[n]])*np.cos(m_2chi[dpoints[n]]))
    z=np.linspace(np.sin(m_2chi[dpoints[n]]),
                  np.sin(m_2chi[dpoints[n]]))
    ax.plot(x,y,z, lw=0.25, color='blue', alpha=1)
    
    x=np.linspace(np.cos(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]]-m_2chi_err[dpoints[n]]),
                  np.cos(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]]+m_2chi_err[dpoints[n]]))
    y=np.linspace(np.sin(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]]-m_2chi_err[dpoints[n]]),
                  np.sin(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]]+m_2chi_err[dpoints[n]]))
    z=np.linspace(np.sin(m_2chi[dpoints[n]]-m_2chi_err[dpoints[n]]),
                  np.sin(m_2chi[dpoints[n]]+m_2chi_err[dpoints[n]]))
    ax.plot(x,y,z, lw=0.25, color='blue', alpha=1)
    
    S3=np.sin(p_2chi[dpoints[n]])
    S2=np.sin(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]])
    S1=np.cos(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]])
    ax.add_artist(Arrow3D([0, S1/np.linalg.norm([S1,S2,S3])],
                          [0, S2/np.linalg.norm([S1,S2,S3])],
                          [0, S3/np.linalg.norm([S1,S2,S3])], mutation_scale=5,
                          lw=0.5, arrowstyle="-|>", color="orange"))
    
    x=np.linspace(np.cos(p_2psi[dpoints[n]]-p_2psi_err[dpoints[n]])*np.cos(p_2chi[dpoints[n]]),
                  np.cos(p_2psi[dpoints[n]]+p_2psi_err[dpoints[n]])*np.cos(p_2chi[dpoints[n]]))
    y=np.linspace(np.sin(p_2psi[dpoints[n]]-p_2psi_err[dpoints[n]])*np.cos(p_2chi[dpoints[n]]),
                  np.sin(p_2psi[dpoints[n]]+p_2psi_err[dpoints[n]])*np.cos(p_2chi[dpoints[n]]))
    z=np.linspace(np.sin(p_2chi[dpoints[n]]),
                  np.sin(p_2chi[dpoints[n]]))
    ax.plot(x,y,z, lw=0.25, color='orange', alpha=1)

    x=np.linspace(np.cos(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]]-p_2chi_err[dpoints[n]]),
                  np.cos(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]]+p_2chi_err[dpoints[n]]))
    y=np.linspace(np.sin(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]]-p_2chi_err[dpoints[n]]),
                  np.sin(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]]+p_2chi_err[dpoints[n]]))
    z=np.linspace(np.sin(p_2chi[dpoints[n]]-p_2chi_err[dpoints[n]]),
                  np.sin(p_2chi[dpoints[n]]+p_2chi_err[dpoints[n]]))
    ax.plot(x,y,z, lw=0.25, color='orange', alpha=1)
    
# Turn off the axis planes
ax.set_axis_off()
plt.show()

#############################################################################
# plotting hemispheres

fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1,2,1, projection='3d')
plot_sphere(ax, arrows='xy', equatorial=False)

ax2 = fig.add_subplot(1,2,2, projection='3d')
plot_sphere(ax2, arrows='xy', equatorial=False)

# Plotting selected datapoints
npoints=200
dpoints=[]
for n in range(npoints):
    dpoints.append(int(random.random()*len(m_dops)))
dpoints=np.array(dpoints)

size = 3
for n in range(npoints):
    S3=np.sin(m_2chi[dpoints[n]])
    S2=np.sin(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]])
    S1=np.cos(m_2psi[dpoints[n]])*np.cos(m_2chi[dpoints[n]])
    if S3 >=0:
        ax.scatter(S1, S2, S3, color='blue', s=size)
    if S3 <=0:
        ax2.scatter(S1, S2, S3, color='blue', s=size)
        
    S3=np.sin(p_2chi[dpoints[n]])
    S2=np.sin(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]])
    S1=np.cos(p_2psi[dpoints[n]])*np.cos(p_2chi[dpoints[n]])
    if S3 >=0:        
        ax.scatter(S1, S2, S3, color='orange', s=size)
    if S3 <=0:        
        ax2.scatter(S1, S2, S3, color='orange', s=size)

#mesh on sphere
for phi in [np.pi/2, (np.pi/2)/3, 2*(np.pi/2)/3]:
    theta = np.linspace(0, 2*np.pi, 100)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    ax.plot(x, y, -z, '--', dashes=(10, 10), color='red', lw=0.5)
    ax2.plot(x, y, z, '--', dashes=(10, 10), color='red', lw=0.5)

for theta in np.linspace(0, 2*np.pi, 12+1):
    phi = np.linspace(np.pi/2, np.pi, 100)
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    ax.plot(x, y, z, '--', dashes=(10, 10), color='red', lw=0.5)
    ax2.plot(x, y, z, '--', dashes=(10, 10), color='red', lw=0.5)

ax.set_axis_off()
ax.view_init(90, 0)
ax.scatter(0,0,1.1,s=10,color='black')
ax.text(0, 0.1, 1.1, '$S_3$', fontweight='bold')
ax2.set_axis_off()
ax2.view_init(-90, 0)
ax2.scatter(0,0,-1.1,marker="x", s=30,color='black')
ax2.text(0, 0.1, 1.1, '$S_3$', fontweight='bold')
plt.show()