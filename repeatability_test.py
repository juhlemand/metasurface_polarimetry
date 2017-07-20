import csv, os
import numpy as np
import matplotlib.pyplot as plt

directory='acquisition\\data\\repeatability' #data location folder
polarimeter_file='polarimeter.txt' #first measurement
os.chdir(directory)

class polarimeter:
    def __init__(self, N_measurements):
        self.N_measurements=N_measurements
        self.data=np.zeros((N_measurements,4))
        self.stdev=np.zeros((N_measurements,4))
        self.poincare=np.zeros((N_measurements,2))
        
    def add(self, index, stokes_vector, stokes_stdev, poincare):
        assert stokes_vector.shape==(4,)
        assert stokes_stdev.shape==(4,)
        assert poincare.shape==(2,)
        self.data[index]=stokes_vector
        self.stdev[index]=stokes_stdev
        self.poincare[index]=poincare

    def __str__(self):
        print(self.data)
        return str(self.stdev)

###################################################################
#parsing set of measurements       
polarimeter_raw=[]
with open(polarimeter_file, 'r') as csvfile:
    reader = csv.reader(csvfile,delimiter=',')
    for row in reader:
        #print(', '.join(row))
        polarimeter_raw.append(row)

N_measurements=0 #number of measurements which have been taken
for row in polarimeter_raw:
    if row[:15]==['#####START#####']:
        N_measurements+=1

polarimeter_data1=polarimeter(N_measurements//2)
polarimeter_data2=polarimeter(N_measurements//2)

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
        s_poincare=np.array([np.mean(np.array(p),0)[3], np.mean(np.array(p),0)[4]])

        if index<N_measurements//2:
            polarimeter_data1.add(index, s_v, s_std, s_poincare)
        elif index>=N_measurements//2:
            polarimeter_data2.add(index-N_measurements//2, s_v, s_std, s_poincare)
        index+=1
        p=[]
        
    elif is_measurement:
        p.append(np.array(row,dtype=float))
 
###############################################################
#%% Analysing and plotting comparison
p_dops1=polarimeter_data1.data.transpose()[0]
p_dops2=polarimeter_data2.data.transpose()[0]

plt.errorbar(range(0,N_measurements//2),p_dops1, alpha=0.5, yerr=polarimeter_data1.stdev.transpose()[0], label='First measurement',fmt='.')
plt.errorbar(range(0,N_measurements//2),p_dops2, alpha=0.5, yerr=polarimeter_data2.stdev.transpose()[0], label='Second measurement',fmt='.')
plt.legend()
plt.show()

f, axarr  = plt.subplots(2,3)
axarr[0][0].scatter(p_dops1, p_dops2,alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][0].plot([0,1],[0,1],alpha=0.3,color='black')
axarr[0][0].set_title('DOP')
axarr[0][0].set_xlabel('First measurement')
axarr[0][0].set_ylabel('Second measurement')

diffs=p_dops1-p_dops2
axarr[1][0].hist(diffs,bins=np.arange(min(diffs), max(diffs) + 0.001, 0.001))
axarr[1][0].axvline(0.0,color='black', alpha=0.25)

#############################################################
#poincare sphere coordinates in radians
p_2psi1=2*polarimeter_data1.poincare.transpose()[0]
p_2psi2=2*polarimeter_data2.poincare.transpose()[0]

p_2chi1=2*polarimeter_data1.poincare.transpose()[1]
p_2chi2=2*polarimeter_data2.poincare.transpose()[1]

axarr[0][1].scatter(p_2psi1, p_2psi2, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][1].set_title('Azimuth $2\psi$')
axarr[0][1].set_xlabel('First measurement (deg)')
axarr[0][1].set_ylabel('Second measurement (deg)')
diffs=p_2psi1-p_2psi2
axarr[1][1].hist(diffs,bins=100)#np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][1].axvline(0.0,color='black', alpha=0.25)

axarr[0][2].scatter(p_2chi1, p_2chi2, alpha=0.5,s=2)#,c=np.arange(0,len(polarimeter_dops)), cmap='viridis')
axarr[0][2].set_title('Altitude $2\chi$')
axarr[0][2].set_xlabel('First measurement (deg)')
axarr[0][2].set_ylabel('Second measurement (deg)')
diffs=p_2chi1-p_2chi2
axarr[1][2].hist(diffs,bins=100)#np.arange(min(diffs), max(diffs) + 0.005, 0.005))
axarr[1][2].axvline(0.0,color='black', alpha=0.25)

plt.show()

