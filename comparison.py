import csv, os
import numpy as np
import matplotlib.pyplot as plt

directory='2017_07_11\\comparison\\'
polarimeter_file='polarimeter.txt'
N_per_polarimeter_meas=100

os.chdir(directory)

def determine_dop(pol_state):
    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]
        
polarimeter_raw=[]
with open(polarimeter_file, 'r') as csvfile:
    reader=csv.reader(csvfile, delimiter=',')
    for row in reader:
        polarimeter_raw.append([])
        for el in row:
            polarimeter_raw[-1].append(float(el))

for i in range(len(polarimeter_raw)//N_per_polarimeter_meas):
    polarimeter_raw[i]=np.average(polarimeter_raw[i*N_per_polarimeter_meas:(i+1)*N_per_polarimeter_meas],0)

polarimeter_raw=np.array(polarimeter_raw[0:len(polarimeter_raw)//N_per_polarimeter_meas])
polarimeter_data=[]
for i in range(len(polarimeter_raw)):
    polarimeter_data.append([polarimeter_raw[i][-2]/N_per_polarimeter_meas, polarimeter_raw[i][0], polarimeter_raw[i][1], polarimeter_raw[i][2]])
polarimeter_data=np.array(polarimeter_data)

fnames=os.listdir()
fnames.remove(polarimeter_file)
fnames=sorted(fnames, key=lambda item: (int(item.partition('_')[0])
                               if item[0].isdigit() else float('inf'), item))
if len(fnames)!=len(polarimeter_data):
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
polarimeter_dops=polarimeter_data.transpose()[0]

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
