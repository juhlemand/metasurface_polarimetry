"""
Created on Wed Jun 21 09:42:59 2017

This is a script to analyze polarimetry calibration data.

@contributors: Noah, Ruoping
"""

import csv, os, pickle
import numpy as np
import matplotlib.pyplot as plt

directory='acquisition/data/calibration4/comparison' #data location folder
polarimeter_file = 'polarimeter.txt' #polarimeter data file
os.chdir(directory)
N_measurements=len(os.listdir())-1 #number of measurements which have been taken

#instrument matrix from calibration
Ainv=np.loadtxt('../Ainv.txt')
Ainv_cov=pickle.load(open( "../Ainv_cov.p", "rb" ))

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
#%% Analyzing and plotting comparison

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

f, axarr  = plt.subplots(2,3)
axarr[0][0].scatter(p_dops, m_dops,alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].errorbar(p_dops, m_dops, xerr=p_dops_err, yerr=m_dops_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].plot([0,1],[0,1],alpha=0.75,color='black')
axarr[0][0].set_title('DOP')
axarr[0][0].set_xlabel('Polarimeter measurement')
axarr[0][0].set_ylabel('Metasurface measurement')
axarr[0][0].set_xlim([-0.1,1.1])
axarr[0][0].set_ylim([-0.1,1.1])

diffs=100*(m_dops-p_dops)  # relative dop error
axarr[1][0].hist(diffs, bins=np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][0].axvline(0.0,color='black', alpha=0.25)
#axarr[1][0].set_title('DOP error metasurface-polarimeter')

##############################################################
#poincare sphere coordinates in radians

# azimuth psi
m_2psi=np.arctan(S2/S1)
p_2psi=np.arctan(polarimeter_data.data.transpose()[2]/polarimeter_data.data.transpose()[1])
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
m_2chi=np.arctan(S3/np.sqrt(S1**2+S2**2))
#https://www.wolframalpha.com/input/?i=derivative+of+atan(x3%2Fsqrt(x1**2%2Bx2**2))
m_2chi_err=np.sqrt((dS1*S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS2*S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS3*np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))**2
                   +2*cov_S2_S1*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S1*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S2*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   )
p_2chi_err=2*np.pi*polarimeter_data.poincare_std.transpose()[2]/180

#taking the first few points as azimuth normalization
az_offset = np.mean(p_2psi[:int(0.05*N_measurements)]-m_2psi[:int(0.05*N_measurements)])

axarr[0][1].scatter(p_2psi, m_2psi, alpha=0.5, s=2.)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].errorbar(p_2psi, m_2psi, xerr=p_2psi_err, yerr=m_2psi_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].plot([np.min(p_2psi), np.max(p_2psi)],[np.min(p_2psi),np.max(p_2psi)],alpha=0.75,color='black')
axarr[0][1].set_title('Azimuth $2\psi$')
axarr[0][1].set_xlabel('Polarimeter measurement (radians)')
axarr[0][1].set_ylabel('Metasurface measurement (radians)')

#diffs=(m_2psi-p_2psi)/(0.5*(m_2psi+p_2psi))
diffs=m_2psi-p_2psi
axarr[1][1].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.01, 0.01))
axarr[1][1].axvline(0.0,color='black', alpha=0.25)

axarr[0][2].scatter(p_2chi, m_2chi, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].errorbar(p_2chi, m_2chi, xerr=p_2chi_err, yerr=m_2chi_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].plot([np.min(p_2chi), np.max(p_2chi)],[np.min(p_2chi),np.max(p_2chi)],alpha=0.75,color='black')
axarr[0][2].set_title('Altitude $2\chi$')
axarr[0][2].set_xlabel('Polarimeter measurement (radians)')
axarr[0][2].set_ylabel('Metasurface measurement (radians)')
#diffs=(m_2chi-p_2chi)/(0.5*(m_2chi+p_2chi))
diffs=m_2chi-p_2chi
axarr[1][2].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.01, 0.01))
axarr[1][2].axvline(0.0,color='black', alpha=0.25)

plt.show()
