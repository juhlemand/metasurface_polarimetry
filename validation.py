import csv, os
import numpy as np
import matplotlib.pyplot as plt

directory='acquisition\\data\\calibration2\\comparison\\' #data location folder
polarimeter_file='polarimeter.txt' #polarimeter data file
N_measurements=100 #number of measurements which have been taken

#instrument matrix from calibration
Ainv=np.array([[ 0.03899718,  0.07547504,  0.04148296,  0.09468381],
               [ 0.12294373, -0.00650341,  0.01732834, -0.19932982],
               [-0.14517754, -0.01293949,  0.12647456, -0.03865285],
               [-0.06027172,  0.20364724, -0.02019657, -0.11410314]])

Ainv_err_rel=np.array([[ 0.11481761,  0.15443711,  0.06899178,  0.29425531],
                       [ 0.0043694 ,  0.00988064,  0.00404934,  0.19303621],
                       [ 0.03765362,  0.01228509,  0.0028296 ,  0.20038666],
                       [ 0.01515934,  0.02031602,  0.015486  ,  0.01874301]])
#https://www.ruhr-uni-bochum.de/ika/forschung/forschungsbereich_kolossa/Daten/Buchkapitel_Uncertainty.pdf

os.chdir(directory)

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

metasurface_data, err_m,=[],[]
err_m2=[]
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
    #covariance matrix-based error calculation, from wikipedia
    cov_m = np.diag(np.std(temp,0)**2)
    #covariance of stokes vector
    cov_stokes = np.dot(np.dot(Ainv,cov_m),np.linalg.inv(Ainv))
    #taking the diagonal as standard error of stokes vector
    err=np.sqrt(np.diagonal(cov_stokes))
    #we also another set of errors from Ainv
    err_m.append(np.sqrt(err**2+np.dot(Ainv_err_rel**2,err**2)))
    #err_m2.append(np.sqrt(np.dot(Ainv_err_rel**2,err**2))
        
err_m=np.array(err_m)
metasurface_data=np.array(metasurface_data)

###############################################################
#%% Analysing and plotting comparison

for i in range(len(metasurface_data)):
    metasurface_data[i]=np.dot(Ainv, metasurface_data[i].transpose())

m_dops=np.sqrt(metasurface_data.transpose()[1]**2+metasurface_data.transpose()[2]**2+metasurface_data.transpose()[3]**2)/metasurface_data.transpose()[0]
p_dops=polarimeter_data.data.transpose()[0]

#error in dop
plt.errorbar(range(0,len(fnames)),m_dops, alpha=0.5, yerr=err_m.transpose()[0], label='metasurface',fmt='.')
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

diffs=m_dops-p_dops
axarr[1][0].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][0].axvline(0.0,color='black', alpha=0.25)
#axarr[1][0].set_title('DOP error metasurface-polarimeter')

##############################################################
#poincare sphere coordinates in radians
p_2psi=np.arctan(polarimeter_data.data.transpose()[2]/polarimeter_data.data.transpose()[1])
m_2psi=np.arctan(metasurface_data.transpose()[2]/metasurface_data.transpose()[1])

p_2chi=np.arctan(polarimeter_data.data.transpose()[3]/np.sqrt(polarimeter_data.data.transpose()[1]**2+polarimeter_data.data.transpose()[2]**2))
m_2chi=np.arctan(metasurface_data.transpose()[3]/np.sqrt(metasurface_data.transpose()[1]**2+metasurface_data.transpose()[2]**2))

#taking the first few points as azimuth normalization
az_offset = np.mean(p_2psi[:int(0.05*N_measurements)]-m_2psi[:int(0.05*N_measurements)])

axarr[0][1].scatter(p_2psi, m_2psi, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].set_title('Azimuth $2\psi$')
axarr[0][1].set_xlabel('Polarimeter measurement (radians)')
axarr[0][1].set_ylabel('Metasurface measurement (radians)')
diffs=m_2psi-p_2psi
axarr[1][1].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][1].axvline(0.0,color='black', alpha=0.25)

axarr[0][2].scatter(-p_2chi, m_2chi, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].set_title('Altitude $2\chi$')
axarr[0][2].set_xlabel('Polarimeter measurement (radians)')
axarr[0][2].set_ylabel('Metasurface measurement (radians)')
diffs=m_2chi+p_2chi
axarr[1][2].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][2].axvline(0.0,color='black', alpha=0.25)


plt.show()
