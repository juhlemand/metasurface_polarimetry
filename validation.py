import csv, os
import numpy as np
import matplotlib.pyplot as plt

directory='acquisition\\data\\comparison\\'
polarimeter_file='polarimeter.txt'
N_measurements=3

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
        
def determine_dop(pol_state):
    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]
        
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
    if row==['#####START#####']:
        is_measurement=1
              
    elif row==['#####END#####']:
        is_measurement=0
        s_v=np.array([np.mean(np.array(p),0)[0],np.mean(np.array(p),0)[1],np.mean(np.array(p),0)[2],np.mean(np.array(p),0)[-2]])
        s_std=np.array([np.std(np.array(p),0)[0],np.std(np.array(p),0)[1],np.std(np.array(p),0)[2],np.std(np.array(p),0)[-2]])
        polarimeter_data.add(index, s_v, s_std)
        index+=1
        p=[]
        
    elif is_measurement:
        p.append(np.array(row,dtype=float))

        
fnames=os.listdir()
fnames.remove(polarimeter_file)
fnames=sorted(fnames, key=lambda item: (int(item.partition('_')[0])
                               if item[0].isdigit() else float('inf'), item))

if len(fnames)!=polarimeter_data.N_measurements:
    raise ValueError

metasurface_data=[]
for n in range(len(fnames)):
    temp=[]
    with open(fnames[n], 'r') as csvfile:
        reader=csv.reader(csvfile, delimiter=',')
        for row in reader:
            temp.append([])
            for el in row:
                temp[-1].append(float(el))
    metasurface_data.append(np.average(temp,0))

metasurface_data=np.array(metasurface_data)

Ainv=np.array([[ 0.0617117 ,  0.07851792,  0.03461167,  0.07703079],
       [ 0.19654736, -0.11104338, -0.07040607,  0.01178421],
       [-0.02460813, -0.15807655,  0.11181989, -0.01704692],
       [ 0.06109838,  0.11286954,  0.02068727, -0.2059895 ]])

for i in range(len(metasurface_data)):
    metasurface_data[i]=np.dot(Ainv, metasurface_data[i].transpose())
    metasurface_data[i][0]=determine_dop(metasurface_data[i])

metasurface_dops=metasurface_data.transpose()[0]
polarimeter_dops=polarimeter_data.data.transpose()[-1]/100.

plt.scatter(range(0,len(fnames)),metasurface_dops, alpha=0.5, label='metasurface')
plt.scatter(range(0,len(fnames)),polarimeter_dops, alpha=0.5, label='thorlabs')
plt.plot((0,len(fnames)), (1.0,1.0), color='gray', alpha=0.5)
plt.legend()
plt.show()

diffs=metasurface_dops-polarimeter_dops
plt.hist(diffs)
plt.show()

plt.scatter(polarimeter_dops, metasurface_dops)
plt.show()
