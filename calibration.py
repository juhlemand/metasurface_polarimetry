# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:42:59 2017
Modified on Thurs Jun 22 2017
Modified on Sat Jun 24 2017

This is a script to analyze polarimetry calibration data.

@author: Noah
"""
import os
import fnmatch
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pyplot as plt
from itertools import compress
matplotlib.use('qt5agg')

linear_pol_extension = 'polarizer_only'  # folder for linear pol data
qwp_R = 'qwp_R'  # folder for qwp at first configuration
qwp_L = 'qwp_L'  # folder for qwp at second configuration
partial_pol = 'partial_pol'  # folder location of partial pol data
comparison = 'polarimeter_comparison'  # folder for comparing polarimeter data

power_meter_error = 0.01

os.chdir('acquisition\\data\\calibration1')

#%% Extract and fit linear polarizer data.

angles = []  # angles of the linear polarizer
inc_powers = []  # incident powers during measurement
pd1_voltage = []  # photodiode 1 voltages
pd2_voltage = []
pd3_voltage = []
pd4_voltage = []
pd1_voltage_err = []  # photodiode errors
pd2_voltage_err = []
pd3_voltage_err = []
pd4_voltage_err = []
fit_functions = []
fit_parameters = []
variances = []

os.chdir(linear_pol_extension)  # go get the linear pol data

for file in os.listdir():
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
            pd1_voltage_err.append(np.std(my_data[:, 0]))
            pd2_voltage_err.append(np.std(my_data[:, 1]))
            pd3_voltage_err.append(np.std(my_data[:, 2]))
            pd4_voltage_err.append(np.std(my_data[:, 3]))
            
        except ValueError: # don't do anything with invalid file name
            pass
       
# rearrange all data as sorted by pol_angle
sorted_lists = sorted(zip(angles, inc_powers, pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage))
# now recover each individual list
angles, inc_powers, pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage =  [[x[i] for x in sorted_lists] for i in range(6)]
        
# begin plotting of obtained data
        
num_angles = len(angles) # number of angular test points
angles = angles[:(num_angles)//2] # cut out part of array from 0 to 170 deg

# split the photodiode powers into two separate lists

pd1_voltage1 = np.array(pd1_voltage[:(num_angles)//2])
pd2_voltage1 = np.array(pd2_voltage[:(num_angles)//2])
pd3_voltage1 = np.array(pd3_voltage[:(num_angles)//2])
pd4_voltage1 = np.array(pd4_voltage[:(num_angles)//2])
pd1_voltage2 = np.array(pd1_voltage[(num_angles)//2:])
pd2_voltage2 = np.array(pd2_voltage[(num_angles)//2:])
pd3_voltage2 = np.array(pd3_voltage[(num_angles)//2:])
pd4_voltage2 = np.array(pd4_voltage[(num_angles)//2:])

pd1_voltage_err1 = np.array(pd1_voltage_err[:(num_angles)//2])
pd2_voltage_err1 = np.array(pd2_voltage_err[:(num_angles)//2])
pd3_voltage_err1 = np.array(pd3_voltage_err[:(num_angles)//2])
pd4_voltage_err1 = np.array(pd4_voltage_err[:(num_angles)//2])
pd1_voltage_err2 = np.array(pd1_voltage_err[(num_angles)//2:])
pd2_voltage_err2 = np.array(pd2_voltage_err[(num_angles)//2:])
pd3_voltage_err2 = np.array(pd3_voltage_err[(num_angles)//2:])
pd4_voltage_err2 = np.array(pd4_voltage_err[(num_angles)//2:])

# slice the incident power list in half as well
inc_powers1 = np.array(inc_powers[:(num_angles)//2])
inc_powers2 = np.array(inc_powers[(num_angles)//2:])

#adding error from power meter, normalizing with incident power
pd1_voltage_err1=np.divide(np.sqrt((pd1_voltage_err1)**2+power_meter_error**2), inc_powers1)
pd2_voltage_err1=np.divide(np.sqrt((pd2_voltage_err1)**2+power_meter_error**2), inc_powers1)
pd3_voltage_err1=np.divide(np.sqrt((pd3_voltage_err1)**2+power_meter_error**2), inc_powers1)
pd4_voltage_err1=np.divide(np.sqrt((pd4_voltage_err1)**2+power_meter_error**2), inc_powers1)
pd1_voltage_err2=np.divide(np.sqrt((pd1_voltage_err2)**2+power_meter_error**2), inc_powers1)
pd2_voltage_err2=np.divide(np.sqrt((pd2_voltage_err2)**2+power_meter_error**2), inc_powers1)
pd3_voltage_err2=np.divide(np.sqrt((pd3_voltage_err2)**2+power_meter_error**2), inc_powers1)
pd4_voltage_err2=np.divide(np.sqrt((pd4_voltage_err2)**2+power_meter_error**2), inc_powers1)

# normalize each by the power incident during measurement
pd1_voltage1 = np.divide(pd1_voltage1, inc_powers1)
pd2_voltage1 = np.divide(pd2_voltage1, inc_powers1)
pd3_voltage1 = np.divide(pd3_voltage1, inc_powers1)
pd4_voltage1 = np.divide(pd4_voltage1, inc_powers1)
pd1_voltage2 = np.divide(pd1_voltage2, inc_powers2)
pd2_voltage2 = np.divide(pd2_voltage2, inc_powers2)
pd3_voltage2 = np.divide(pd3_voltage2, inc_powers2)
pd4_voltage2 = np.divide(pd4_voltage2, inc_powers2)

# now average the two
pd1_voltage = (pd1_voltage1+pd1_voltage2)/2
pd2_voltage = (pd2_voltage1+pd2_voltage2)/2
pd3_voltage = (pd3_voltage1+pd3_voltage2)/2
pd4_voltage = (pd4_voltage1+pd4_voltage2)/2
pd1_voltage_err = np.sqrt((pd1_voltage_err1**2+pd1_voltage_err2**2))/2
pd2_voltage_err = np.sqrt((pd2_voltage_err1**2+pd2_voltage_err2**2))/2
pd3_voltage_err = np.sqrt((pd3_voltage_err1**2+pd3_voltage_err2**2))/2
pd4_voltage_err = np.sqrt((pd4_voltage_err1**2+pd4_voltage_err2**2))/2

# construct a larger matrix to hold the voltages
pd_voltages = np.vstack((pd1_voltage, pd2_voltage, pd3_voltage, pd4_voltage))
pd_errs = np.vstack((pd1_voltage_err,pd2_voltage_err,pd3_voltage_err,pd4_voltage_err))
angles = np.array(angles)
inc_powers = np.array(inc_powers)

# define a fitting function
def fit_function(theta, a, b, c):
    return a + b*np.cos(2*theta*np.pi/180) + c*np.sin(2*theta*np.pi/180)

plt.figure
plt.yticks([]) # y units are arbitrary

thetas = np.linspace(0, 180, 1000)

for i in range(0,4):
    x = angles
    y = pd_voltages[i, :]
    err = pd_errs[i, :]
    popt, variance = curve_fit(fit_function, x, y)
    
    #standard error
    standard_err = np.sqrt(np.diag(variance))
    print(standard_err)
    
    plt.plot(thetas, fit_function(thetas, *popt), linewidth=0.5)
    plt.errorbar(x, y, fmt=" ", yerr=err)
    
    # store the fits as anonymous functions
    fit_functions.append(lambda theta: fit_function(theta, *popt))
    fit_parameters.append(popt)
    variances.append(variance)
plt.show()

#%% move onto the qwpR part of the calibration

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
    if i == 1:
        file = qwp_R
    else:
        file = qwp_L
    
    os.chdir(file)  # go get the qwp+pol data
    for file in os.listdir():
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
                pass
    
    pol_anglesR = np.array(pol_anglesR)
    qwp_anglesR = np.array(qwp_anglesR)
    qwp_power_incR = np.array(qwp_power_incR)
    
    numAngles = len(pol_anglesR)
    
    # normalize the photodiode voltages
    pd1_voltageQR = np.divide(np.array(pd1_voltageQR), qwp_power_incR)
    pd2_voltageQR = np.divide(np.array(pd2_voltageQR), qwp_power_incR)
    pd3_voltageQR = np.divide(np.array(pd3_voltageQR), qwp_power_incR)
    pd4_voltageQR = np.divide(np.array(pd4_voltageQR), qwp_power_incR)
    
    # average all the values in the list
    
    if i == 1:
        pd1R = np.mean(pd1_voltageQR)
        pd2R = np.mean(pd2_voltageQR)
        pd3R = np.mean(pd3_voltageQR)
        pd4R = np.mean(pd4_voltageQR)
    else:
        pd1L = np.mean(pd1_voltageQR)
        pd2L = np.mean(pd2_voltageQR)
        pd3L = np.mean(pd3_voltageQR)
        pd4L = np.mean(pd4_voltageQR)

#%% Construct the instrument matrix
A3 = np.array([(pd1R-pd1L)/2, (pd2R-pd2L)/2, (pd3R-pd3L)/2, (pd4R-pd4L)/2])
A3 = np.matrix.transpose(np.array([A3]))
A02 = np.array(fit_parameters) # data from linear calibration
A = np.hstack((A02, A3)) # This is the instrument matrix!

Ainv = np.linalg.inv(A) # this is the inverse of the instrument matrix
print('Instrument matrix Ainv:')
print(Ainv)

A_cond=np.linalg.cond(A, p=2)#condition number
print('Condition number:')
print(A_cond)

# now let's define a function to reonstruct polarization state
# measurement is a four vector of measured intensities
def determine_stokes(measurement):
    try:
        return np.dot(Ainv, measurement)
    except ValueError:
        raise('Input is not a 4 x 1 intensity measurement!')

def determine_dop(pol_state):
    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]
        
#%% Now check to see that results of instrument matrix make sense

dops = []
stokes_list = []
for i in range(len(pd1_voltage)):
    stokes = np.dot(Ainv, pd_voltages[:, i])
    dop = np.sqrt(stokes[1]**2 + stokes[2]**2 + stokes[3]**2)/stokes[0]
    dops.append(dop)
    stokes_list.append(stokes/stokes[0])
    
plt.figure(2)
plt.scatter(np.arange(len(dops)),dops)
plt.plot([0,len(dops)],[1,1],alpha=0.3, color='black')
plt.show()

L = np.array([pd1L, pd2L, pd3L, pd4L])
stokes=np.dot(Ainv, L)

#%% Extract and analyze the partial pol data

pol_angles2=[]
pd1_partialV = []
pd2_partialV = []
pd3_partialV = []
pd4_partialV = []

os.chdir('..')  # move up a level
os.chdir(partial_pol)  # go get the qwp+pol data
for file in os.listdir():
    if fnmatch.fnmatch(file, '*.txt'):
        params = file.split('_')
        try:  # if file name matches format do stuff with it
            pol_angle = float(params[1][:-7])  # get rid of 'deg' in string
            pol_angles2.append(pol_angle)
                
            # make file into array
            data = np.genfromtxt(file, delimiter=',')
            pd1_partialV.append(np.mean(data[:, 0]))
            pd2_partialV.append(np.mean(data[:, 1]))
            pd3_partialV.append(np.mean(data[:, 2]))
            pd4_partialV.append(np.mean(data[:, 3]))
        except ValueError:  # don't do anything with invalid file name
            pass
        
# rearrange all data as sorted by pol_angles2
sorted_lists = sorted(zip(pol_angles2, pd1_partialV, pd2_partialV, pd3_partialV, pd4_partialV))

# now recover each individual list
pol_angles2, pd1_partialV, pd2_partialV, pd3_partialV, pd4_partialV =  [[x[i] for x in sorted_lists] for i in range(5)]

num_angles = len(pol_angles2)
pol_angles2 = pol_angles2[:num_angles//2]

pd1_partialV = np.divide(pd1_partialV[:num_angles//2] + pd1_partialV[num_angles//2:], 2)
pd2_partialV = np.divide(pd2_partialV[:num_angles//2] + pd2_partialV[num_angles//2:], 2)
pd3_partialV = np.divide(pd3_partialV[:num_angles//2] + pd3_partialV[num_angles//2:], 2)
pd4_partialV = np.divide(pd4_partialV[:num_angles//2] + pd4_partialV[num_angles//2:], 2)
            
partial_dops = np.zeros(len(pol_angles2))
for i in range(len(pol_angles2)):
    i_measured = np.array([pd1_partialV[i], pd2_partialV[i], pd3_partialV[i], pd4_partialV[i]])
    stokes_temp = np.dot(Ainv, i_measured)
    
    partial_dops[i] = np.sqrt(stokes_temp[1]**2 + stokes_temp[2]**2 + stokes_temp[3]**2)/stokes_temp[0]

plt.figure(3)
plt.plot(pol_angles2, partial_dops, ".", markersize=5)
plt.plot([0,180],[1,1],color='black',alpha=0.25)
plt.xlabel('$\Theta_{LP} (\circ)$', fontsize='12', fontname='Sans Serif')
plt.ylabel('Degree of Polarization (DOP)', fontsize='12')
partial_pol_fig = plt.gca()
partial_pol_fig.tick_params(axis='x', labelsize=16, direction='out', length=5)
partial_pol_fig.set_ylim([0, 1.05])


thetas = np.linspace(0, np.max(pol_angles2), 1000)
def partial_pol_func(x, scale, offset):
    return scale * np.abs(np.cos(2*(x-offset)*np.pi/180))

popt2, variance2 = curve_fit(partial_pol_func, pol_angles2, partial_dops)
plt.plot(thetas, partial_pol_func(thetas, *popt2), linewidth = 1, alpha=0.5)
plt.show()
# now save the figure
#file_name = os.getcwd() + 'partial_pol.svg'
#plt.savefig(file_name, format='svg')

##
###%% Analyze data vs. the Thorlabs polarimeter.
##
##os.chdir('..')  # move up a level
##os.chdir(comparison)  # cd into the folder with comparison measurements
##files = os.listdir()
##txp_files = [fnmatch.fnmatch(file, 'TXP_*.csv') for file in files]
##txp_files = compress(files, txp_files)
##measured_files = [fnmatch.fnmatch(file, 'pol_*.txt') for file in files]
##measured_files = compress(files, measured_files)
##txp_files = [f for f in txp_files] # exhaust the iterators
##measured_files = [f for f in measured_files]
### now sort the list by measurement number
##txp_files = sorted(txp_files, key=lambda x: int((x.split('_')[1][:-4])))
##measured_files = sorted(measured_files, key=lambda x: int((x.split('_')[1][:-4])))
##
##pd_voltages = np.zeros((len(measured_files), 4))  # a place to store the measurements
##for file, i in zip(measured_files, range(len(measured_files))):
##    data = np.genfromtxt(file, delimiter=',')
##    # append the data values to the array
##    data = np.mean(data, axis=0)
##    pd_voltages[i, :] = np.array([data[0], data[1], data[2], data[3]])
##
### loop over files from the polarimeter
##s1s = np.zeros(len(txp_files))
##s2s = np.zeros(len(txp_files))
##s3s = np.zeros(len(txp_files))
##azs = np.zeros(len(txp_files))
##ells = np.zeros(len(txp_files))
##dops = np.zeros(len(txp_files))
##pwrs = np.zeros(len(txp_files))
##for file, i in zip(txp_files, range(len(txp_files))):
##    with open(file, mode='rt') as csvfile:
##        # try to fetch data from the file but even if error assure it's closed        
##        try:        
##            # go through the csv header which contains no information     
##            for ii in range(23):
##                csvfile.__next__()
##    
##            s1sTemp = []
##            s2sTemp = []
##            s3sTemp = []
##            azsTemp = []
##            ellsTemp = []
##            dopsTemp = []
##            pwrsTemp = []      
##    
##            # go through the file line-by-line and extract data
##            for data_line in csvfile:
##                recorded = data_line.split(',')
##                s1 = float(recorded[1])
##                s2 = float(recorded[2])
##                s3 = float(recorded[3])
##                az = float(recorded[4])
##                ell = float(recorded[5])
##                dop = float(recorded[8])
##                pwr = float(recorded[10])
##
##    
##                s1sTemp.append(s1)
##                s2sTemp.append(s2)
##                s3sTemp.append(s3)
##                azsTemp.append(az)
##                ellsTemp.append(ell)
##                dopsTemp.append(dop)
##                pwrsTemp.append(pwr)
##            # average the data from the file and add to array of data for all files
##            s1s[i] = np.mean(np.array(s1sTemp))
##            s2s[i] = np.mean(np.array(s2sTemp))
##            s3s[i] = np.mean(np.array(s3sTemp))
##            azs[i] = np.mean(np.array(azsTemp))
##            ells[i] = np.mean(np.array(ellsTemp))
##            dops[i] = np.mean(np.array(dopsTemp))
##            pwrs[i] = np.mean(np.array(pwrsTemp))
##        finally:
##            csvfile.close()
##
### need to divide dops by 100 from percent form
##dops = dops/100
##
### analyze the azimuth reference measurement
##reference_file = 'pol_reference_TXP.csv'
##reference_file_reader = open(file, mode='rt')
##reference_az = []
### fetch the data from the file but make sure to close it in case of error
##try:
##    for i in range(23):
##        reference_file_reader.__next__()
##    for data_line in reference_file_reader:
##        recorded = data_line.split(',')
##        reference_az.append(float(recorded[4]))
##    reference_az = np.mean(reference_az)
##finally:
##    reference_file_reader.close()
##
### loop over files from our polarimeter
##pol_states = np.zeros((len(measured_files), 4))
##dop_states = np.zeros(len(measured_files))
##azimuths = np.zeros((len(measured_files)))
##phase = np.zeros(len(measured_files))
##for file, i in zip(measured_files, range(len(measured_files))):            
##    # make file into array
##    my_data = np.genfromtxt(file, delimiter=',')
##    det1 = np.mean(my_data[:, 0])
##    det2 = np.mean(my_data[:, 1])
##    det3 = np.mean(my_data[:, 2])
##    det4 = np.mean(my_data[:, 3])
##    
##    # compute the stokes vector
##    pol_state = determine_stokes(np.array([det1, det2, det3, det4]))
##    # normalize by the first element
##    pol_state = pol_state/pol_state[0]
##    # add to the running vector of stokes vectors
##    pol_states[i, :] = pol_state
##    dop_states[i] = determine_dop(pol_state)
##    azimuths[i] = np.degrees(np.arctan2(pol_state[2], pol_state[1]))
##    equator_component = np.sqrt(pol_state[1]**2+pol_state[2]**2)
##    phase[i] = np.degrees(np.arctan2(pol_state[3], equator_component))
##
##
###%% Plot the polarimeter comparison data:
##plt.figure(4)
### plot data from the ThorLabs polarimeter
##plt.scatter(range(len(dops)), dops, marker='o')
##helicities = pol_states[:, 3]
##plt.show()
### plot data from our polarimeter
##plt.scatter(range(len(dop_states)), dop_states, marker='x')
##plt.show()
