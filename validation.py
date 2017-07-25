"""
Created on Wed Jun 21 09:42:59 2017

This is a script to analyze polarimetry calibration data.

@contributors: Noah, Ruoping
"""

import csv, os, pickle
import numpy as np
import matplotlib.pyplot as plt

directory='acquisition/data/calibration1/comparison1.2' #data location folder
polarimeter_file='polarimeter.txt' #polarimeter data file
os.chdir(directory)
N_measurements=len(os.listdir())-1 #number of measurements which have been taken

#instrument matrix from calibration
Ainv=np.loadtxt('../Ainv.txt')
Ainv_cov=pickle.load(open( "../Ainv_cov.p", "rb" ))

#https://www.ruhr-uni-bochum.de/ika/forschung/forschungsbereich_kolossa/Daten/Buchkapitel_Uncertainty.pdf

def covS(i,j, D, I, Dcov, Icov):
    ''' This function returns the covariance matrix of the result of Ainv*I
    '''
    assert len(I)==4
    assert D.shape==(4,4)
    assert Dcov.shape==(4,4,4,4)
    s=0.0
    for a in range(4):
        for b in range(4):            
            s+=I[a]*I[b]*Dcov[i][a][j][b]
    for k in range(4):
        for l in range(4):
            s+=D[i][k]*D[j][l]*Icov[k][l]
    return s

class polarimeter:
    def __init__(self, N_measurements):
        self.N_measurements=N_measurements
        self.data=np.zeros((N_measurements,4))
        self.stdev=np.zeros((N_measurements,4))
        
    def add(self, index, stokes_vector, stokes_stdev):
        assert stokes_vector.shape==(4,)
        assert stokes_stdev.shape==(4,)
        self.data[index]=stokes_vector
        self.stdev[index]=stokes_stdev

    def __str__(self):
        print(self.data)
        return str(self.stdev)
        
polarimeter_raw=[]
with open(polarimeter_file, 'r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    for row in reader:
        #print(', '.join(row))
        polarimeter_raw.append(row)
polarimeter_raw=np.array(polarimeter_raw)

#####################################################
#%% Parsing polarimeter file
polarimeter_data=polarimeter(N_measurements)
p=[]
is_measurement=0
index=0
for row in polarimeter_raw:
    if row[:15]==['#####START#####']:
        is_measurement=1
              
    elif row[:13]==['#####END#####']:
        is_measurement=0
        s_v=np.array([0.01*np.mean(np.array(p),0)[-2], np.mean(np.array(p),0)[0],np.mean(np.array(p),0)[1],np.mean(np.array(p),0)[2]])
        s_std=0.01*np.array([np.std(np.array(p),0)[-2], np.std(np.array(p),0)[0],np.std(np.array(p),0)[1],np.std(np.array(p),0)[2]])
        polarimeter_data.add(index, s_v, s_std)
        index+=1
        p=[]
        
    elif is_measurement:
        p.append(np.array(row,dtype=float))

###############################################################
#%% Parsing metasurface data
fnames=os.listdir()
fnames.remove(polarimeter_file)
fnames=sorted(fnames, key=lambda item: (int(item.partition('_')[0])
                               if item[0].isdigit() else float('inf'), item))

if len(fnames)!=polarimeter_data.N_measurements:
    raise ValueError

metasurface_data, cov_m,=[],[]
cov_m=[]
for n in range(len(fnames)):
    temp=[]
    with open(fnames[n], 'r') as csvfile:
        reader=csv.reader(csvfile, delimiter=',')
        for row in reader:
            temp.append([])
            for el in row:
                temp[-1].append(float(el))
    temp=np.array(temp)
    metasurface_data.append(np.average(temp,0))
    I=np.average(temp,0)
    # assume covariance matrix is diagonal
    Icov = np.diag(np.std(temp,0)**2)
    #covariance of stokes vector
    cov_stokes = np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            cov_stokes[i][j]=covS(i,j, Ainv, I, Ainv_cov, Icov)
    #taking the diagonal as standard error of stokes vector
    cov_m.append(cov_stokes)
        
cov_m=np.array(cov_m)
metasurface_data=np.array(metasurface_data)

###############################################################
#%% Analysing and plotting comparison

for i in range(len(metasurface_data)):
    metasurface_data[i]=np.dot(Ainv, metasurface_data[i].transpose())

m_dops=np.sqrt(metasurface_data.transpose()[1]**2+metasurface_data.transpose()[2]**2+metasurface_data.transpose()[3]**2)/metasurface_data.transpose()[0]
p_dops=polarimeter_data.data.transpose()[0]

plt.errorbar(range(0,len(fnames)),m_dops, alpha=0.5, label='metasurface',fmt='.')
plt.errorbar(range(0,len(fnames)),p_dops, alpha=0.5, yerr=polarimeter_data.stdev.transpose()[0], label='thorlabs',fmt='.')
plt.plot((0,len(fnames)), (1.0,1.0), color='gray', alpha=0.5)
plt.ylim([0,1.1])
plt.legend()
plt.show()

f, axarr  = plt.subplots(2,3)
axarr[0][0].scatter(p_dops, m_dops,alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].plot([0,1],[0,1],alpha=0.3,color='black')
axarr[0][0].set_title('DOP')
axarr[0][0].set_xlabel('Polarimeter measurement')
axarr[0][0].set_ylabel('Metasurface measurement')

diffs=(m_dops-p_dops)/(0.5*(m_dops+p_dops))
axarr[1][0].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][0].axvline(0.0,color='black', alpha=0.25)
#axarr[1][0].set_title('DOP error metasurface-polarimeter')

##############################################################
#poincare sphere coordinates in radians

S3=metasurface_data.transpose()[3]
S2=metasurface_data.transpose()[2]
S1=metasurface_data.transpose()[1]
S0=metasurface_data.transpose()[0]
dS3=np.sqrt(cov_m[:,3,3])
dS2=np.sqrt(cov_m[:,2,2])
dS1=np.sqrt(cov_m[:,1,1])
dS0=np.sqrt(cov_m[:,0,0])
cov_S2_S1=cov_m[:,2,1]
cov_S3_S1=cov_m[:,3,1]
cov_S3_S2=cov_m[:,3,2]

#error in dop
# dop = sqrt(S1**2+S2**2+S3**2)/S0
#https://www.wolframalpha.com/input/?i=derivative+of+sqrt(x1**2%2Bx2**2%2Bx3**2)%2Fx0


# azimuth psi
m_2psi=np.arctan(S2/S1)
p_2psi=np.arctan(polarimeter_data.data.transpose()[2]/polarimeter_data.data.transpose()[1])
#p_2psi_err=np.
#stdev error in S2/S1:
#http://123.physics.ucdavis.edu/week_9_files/taylor_209-226.pdf
m_2psi_err=np.sqrt((S1*dS2/(S1**2+S2**2))**2+(S2*dS1/(S1**2+S2**2))**2
                   -2*S1*S2*cov_S2_S1/(S1**2+S2**2)**2)

# altitude chi
p_2chi=np.arctan(polarimeter_data.data.transpose()[3]/np.sqrt(polarimeter_data.data.transpose()[1]**2+polarimeter_data.data.transpose()[2]**2))
m_2chi=np.arctan(S3/np.sqrt(S1**2+S2**2))
#https://www.wolframalpha.com/input/?i=derivative+of+atan(x3%2Fsqrt(x1**2%2Bx2**2))
m_2chi_err=np.sqrt((dS1*S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS2*S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))**2
                   +(dS3*np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))**2
                   +2*cov_S2_S1*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S1*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S1*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   +2*cov_S3_S2*(np.sqrt(S1**2+S2**2)/(S1**2+S2**2+S3**2))*(-S2*S3/(np.sqrt(S1**2+S2**2)*(S1**2+S2**2+S3**2)))
                   )


#taking the first few points as azimuth normalization
az_offset = np.mean(p_2psi[:int(0.05*N_measurements)]-m_2psi[:int(0.05*N_measurements)])

axarr[0][1].scatter(p_2psi, m_2psi, alpha=0.5, s=2.)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].errorbar(p_2psi, m_2psi, yerr=m_2psi_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].set_title('Azimuth $2\psi$')
axarr[0][1].set_xlabel('Polarimeter measurement (radians)')
axarr[0][1].set_ylabel('Metasurface measurement (radians)')

#diffs=(m_2psi-p_2psi)/(0.5*(m_2psi+p_2psi))
diffs=m_2psi-p_2psi
axarr[1][1].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.01, 0.01))
axarr[1][1].axvline(0.0,color='black', alpha=0.25)

axarr[0][2].scatter(p_2chi, m_2chi, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].errorbar(p_2chi, m_2chi, yerr=m_2chi_err, alpha=0.5, fmt=' ')#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].set_title('Altitude $2\chi$')
axarr[0][2].set_xlabel('Polarimeter measurement (radians)')
axarr[0][2].set_ylabel('Metasurface measurement (radians)')
#diffs=(m_2chi-p_2chi)/(0.5*(m_2chi+p_2chi))
diffs=m_2chi-p_2chi
axarr[1][2].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.01, 0.01))
axarr[1][2].axvline(0.0,color='black', alpha=0.25)

plt.show()
