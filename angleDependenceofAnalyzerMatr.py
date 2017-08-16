# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 09:42:59 2017

This is a script to analyze polarimetry calibration data.

@contributors: Noah, Ruoping
"""
#%reset -f

import os, re, pickle
import fnmatch
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
from matplotlib.ticker import AutoMinorLocator

plt.close('all')

linear_pol_extension = 'polarizer_only'  # folder for linear pol data
qwp_R = 'qwp_R'  # folder for qwp at first configuration
qwp_L = 'qwp_L'  # folder for qwp at second configuration

partial_pol = 'partial_pol8'  # folder location of partial pol data

power_meter_error = 0.001 #Error in power meter reading from ambient light, unit in mW


if 'linux' or 'darwin' in sys.platform:
    data_dir = 'acquisition/data/incident_angles_calibration'
else:
    data_dir = 'acquisition\\data\\incident_angles_calibration\\20deg'

os.chdir(data_dir)
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


#%% Extract and fit linear polarizer data.
q=0
A=np.zeros((len(angledirs),4,4))
Ainv=np.zeros((len(angledirs),4,4))
for folder in angledirs:
    os.chdir(folder)

    angles = []  # angles off the linear polarizer
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
                #Here error is the standard deviation of the individual measurements
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
    
    # normalize each by the power incident during measurement
    pd1_voltage1 = pd1_voltage1/inc_powers1
    pd2_voltage1 = pd2_voltage1/inc_powers1
    pd3_voltage1 = pd3_voltage1/inc_powers1
    pd4_voltage1 = pd4_voltage1/inc_powers1
    pd1_voltage2 = pd1_voltage2/inc_powers2
    pd2_voltage2 = pd2_voltage2/inc_powers2
    pd3_voltage2 = pd3_voltage2/inc_powers2
    pd4_voltage2 = pd4_voltage2/inc_powers2
    
    # adding error from power meter, propagating error through normalization
    pd1_voltage_err1 = pd1_voltage1 * np.sqrt((pd1_voltage_err1/pd1_voltage1)**2 + (power_meter_error/inc_powers1)**2)
    pd2_voltage_err1 = pd2_voltage1 * np.sqrt((pd2_voltage_err1/pd2_voltage1)**2 + (power_meter_error/inc_powers1)**2)
    pd3_voltage_err1 = pd3_voltage1 * np.sqrt((pd3_voltage_err1/pd3_voltage1)**2 + (power_meter_error/inc_powers1)**2)
    pd4_voltage_err1 = pd1_voltage1 * np.sqrt((pd4_voltage_err1/pd4_voltage1)**2 + (power_meter_error/inc_powers1)**2)
    pd1_voltage_err2 = pd1_voltage2 * np.sqrt((pd1_voltage_err2/pd1_voltage2)**2 + (power_meter_error/inc_powers2)**2)
    pd2_voltage_err2 = pd2_voltage2 * np.sqrt((pd2_voltage_err2/pd2_voltage2)**2 + (power_meter_error/inc_powers2)**2)
    pd3_voltage_err2 = pd3_voltage2 * np.sqrt((pd3_voltage_err2/pd3_voltage2)**2 + (power_meter_error/inc_powers2)**2)
    pd4_voltage_err2 = pd4_voltage2 * np.sqrt((pd4_voltage_err2/pd4_voltage2)**2 + (power_meter_error/inc_powers2)**2)
    
    
    #pd1_voltage_err1=np.sqrt((pd1_voltage_err1/inc_powers1)**2+(power_meter_error*pd1_voltage1/(inc_powers1*inc_powers1))**2)
    #pd2_voltage_err1=np.sqrt((pd2_voltage_err1/inc_powers1)**2+(power_meter_error*pd2_voltage1/(inc_powers1*inc_powers1))**2)
    #pd3_voltage_err1=np.sqrt((pd3_voltage_err1/inc_powers1)**2+(power_meter_error*pd3_voltage1/(inc_powers1*inc_powers1))**2)
    #pd4_voltage_err1=np.sqrt((pd4_voltage_err1/inc_powers1)**2+(power_meter_error*pd4_voltage1/(inc_powers1*inc_powers1))**2)
    #pd1_voltage_err2=np.sqrt((pd1_voltage_err2/inc_powers1)**2+(power_meter_error*pd1_voltage2/(inc_powers1*inc_powers1))**2)
    #pd2_voltage_err2=np.sqrt((pd2_voltage_err2/inc_powers1)**2+(power_meter_error*pd2_voltage2/(inc_powers1*inc_powers1))**2)
    #pd3_voltage_err2=np.sqrt((pd3_voltage_err2/inc_powers1)**2+(power_meter_error*pd3_voltage2/(inc_powers1*inc_powers1))**2)
    #pd4_voltage_err2=np.sqrt((pd4_voltage_err2/inc_powers1)**2+(power_meter_error*pd4_voltage2/(inc_powers1*inc_powers1))**2)
    
    # now average the two and propagate the error through the averaging
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
    # convert these lists to proper numpy arrays
    angles = np.array(angles)
    inc_powers = np.array(inc_powers)
    
    
    # Now plot linear polarization cal data
    
    # define a function to plot the resulting graph
    def linear_cal_fig(axes, yerror, xdata, ydata, min_angle, max_angle):
        
        # define a fitting function
        def fit_function(theta, a, b, c):
            return a + b*np.cos(2*theta*np.pi/180) + c*np.sin(2*theta*np.pi/180)
        
        axes.set_yticks([]) # y units are arbitrary, so no ticks
        
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
        yerror = yerror[:, mask]
        thetas = np.linspace(min_angle, max_angle, 1000)        
        colors = [(0, 0, 1, 0.5), (0, 0, 0, 0.5), (0, 0.5, 0.25, 0.5), (0.5, 0, 0.5, 0.5)]    
            
        for i in range(0, 4):
            err = np.transpose(yerror[i, :])
            y = np.transpose(ydata[i, :])
            axes.plot(thetas, fit_function(thetas, *fit_parameters[i]), linewidth=1.5, color = colors[i], label='PD #' + str(i+1))
            axes.errorbar(xdata, y, fmt=" ", yerr=err, markersize = 6, ecolor = 'r', color = 'r', capsize=2.5, elinewidth=2.5) 
 #           print(fit_parameters[i])
            
        axes.legend(prop={'size': 9})
        axes.set_ylim([0, 1.02*np.max(pd_voltages)])
        axes.set_xlim([min_angle, max_angle])
        minor_locatorx = AutoMinorLocator(2)
        minor_locatory = AutoMinorLocator(2)
    
        axes.xaxis.set_minor_locator(minor_locatorx)
        axes.yaxis.set_minor_locator(minor_locatory)
    
        major_length = 7
        major_width = 1.5
        minor_length = 5
        minor_width = 1.5    
    
        axes.tick_params(axis='x', labelsize = 14, direction='in', length = major_length, width=major_width, which = 'major', top='off', color='k')
        axes.tick_params(axis='x', labelsize = 14, direction='in', length = minor_length, width=minor_width, which = 'minor', top='off', color='k')
        axes.tick_params(axis='y', labelsize = 14, direction='in', length = major_length, width=major_width, which = 'major', top='off', color='k')
        axes.tick_params(axis='y', labelsize = 14, direction='in', length = minor_length, width=minor_width, which = 'minor', top='off', color='k')
        axes.set_xlim([min_angle, max_angle])
        return fit_errs
    
    
    save_fig = 0
    
    # now plot it
    min_angle = 0
    max_angle = 180
    fig = plt.figure()
    ax = plt.gca()
    fit_errs = linear_cal_fig(ax, pd_errs, angles, pd_voltages, min_angle, max_angle)
    
    if save_fig:
        file_name = 'linear_cal.svg'
        os.chdir('../../../../Graphics')
        plt.savefig(file_name, format='svg')
        os.chdir('..\\' + data_dir + '\\' + linear_pol_extension)
    
    plt.show()
    
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
        pd1_voltage_err = []
        pd2_voltage_err = []
        pd3_voltage_err = []
        pd4_voltage_err = []
        
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
                    pd1_voltage_err.append(np.std(my_data[:, 0]))
                    pd2_voltage_err.append(np.std(my_data[:, 1]))
                    pd3_voltage_err.append(np.std(my_data[:, 2]))
                    pd4_voltage_err.append(np.std(my_data[:, 3]))
                except ValueError:  # don't do anything with invalid file name
                    print('file skipped ', file)
                    pass
        
        pol_anglesR = np.array(pol_anglesR)
        qwp_anglesR = np.array(qwp_anglesR)
        qwp_power_incR = np.array(qwp_power_incR)
        
        numAngles = len(pol_anglesR)
    
        # converting to array
        pd1_voltage_err=np.array(pd1_voltage_err)
        pd2_voltage_err=np.array(pd2_voltage_err)
        pd3_voltage_err=np.array(pd3_voltage_err)
        pd4_voltage_err=np.array(pd4_voltage_err)
        
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
        
        # add in voltage from the power meter
        pd1_voltage_err = pd1_voltageQR * np.sqrt((pd1_voltage_err/pd1_voltageQR)**2 + (power_meter_error/qwp_power_incR)**2)
        pd2_voltage_err = pd2_voltageQR * np.sqrt((pd2_voltage_err/pd2_voltageQR)**2 + (power_meter_error/qwp_power_incR)**2)
        pd3_voltage_err = pd3_voltageQR * np.sqrt((pd3_voltage_err/pd3_voltageQR)**2 + (power_meter_error/qwp_power_incR)**2)
        pd4_voltage_err = pd4_voltageQR * np.sqrt((pd4_voltage_err/pd4_voltageQR)**2 + (power_meter_error/qwp_power_incR)**2)
        
        # adding in the error from power meter
    #    pd1_voltage_err=np.sqrt((pd1_voltage_err/qwp_power_incR)**2+(power_meter_error*pd1_voltage_err/(qwp_power_incR**2))**2)
    #    pd2_voltage_err=np.sqrt((pd2_voltage_err/qwp_power_incR)**2+(power_meter_error*pd2_voltage_err/(qwp_power_incR**2))**2)
    #    pd3_voltage_err=np.sqrt((pd3_voltage_err/qwp_power_incR)**2+(power_meter_error*pd3_voltage_err/(qwp_power_incR**2))**2)
    #    pd4_voltage_err=np.sqrt((pd4_voltage_err/qwp_power_incR)**2+(power_meter_error*pd4_voltage_err/(qwp_power_incR**2))**2)
    
      # average all the values in the list    
        if i == 0:
            # compute the mean value measured on each photodiode
            pd1R = np.mean(pd1_voltageQR)
            pd2R = np.mean(pd2_voltageQR)
            pd3R = np.mean(pd3_voltageQR)
            pd4R = np.mean(pd4_voltageQR)
            # measurement error is standard deviation of sets of 90 deg apart measurements
            pd1R_err = np.std(qwp_err(pd1_voltageQR))
            pd2R_err = np.std(qwp_err(pd2_voltageQR))
            pd3R_err = np.std(qwp_err(pd3_voltageQR))
            pd4R_err = np.std(qwp_err(pd4_voltageQR))
            
            plt.errorbar(pol_anglesR, pd1_voltageQR, yerr=pd1_voltage_err, fmt=' ', color='red')
            plt.errorbar(pol_anglesR,pd2_voltageQR, yerr=pd2_voltage_err, fmt=' ', color='blue')
            plt.errorbar(pol_anglesR, pd3_voltageQR, yerr=pd3_voltage_err, fmt=' ', color='green')
            plt.errorbar(pol_anglesR, pd4_voltageQR,yerr=pd4_voltage_err, fmt=' ', color='orange')
            #print(len(pd1_voltageQR))
        elif i == 1:
            pd1L = np.mean(pd1_voltageQR)
            pd2L = np.mean(pd2_voltageQR)
            pd3L = np.mean(pd3_voltageQR)
            pd4L = np.mean(pd4_voltageQR)
            pd1L_err = np.std(qwp_err(pd1_voltageQR))
            pd2L_err = np.std(qwp_err(pd2_voltageQR))
            pd3L_err = np.std(qwp_err(pd3_voltageQR))
            pd4L_err = np.std(qwp_err(pd4_voltageQR))
            plt.errorbar(pol_anglesR, pd1_voltageQR, yerr=pd1_voltage_err, fmt=' ', color='red', alpha=0.5)
            plt.errorbar(pol_anglesR, pd2_voltageQR, yerr=pd2_voltage_err, fmt=' ', color='blue', alpha=0.5)
            plt.errorbar(pol_anglesR, pd3_voltageQR, yerr=pd3_voltage_err, fmt=' ', color='green', alpha=0.5)
            plt.errorbar(pol_anglesR, pd4_voltageQR, yerr=pd4_voltage_err, fmt=' ', color='orange', alpha=0.5)
            #print(len(pd1_voltageQR))
    plt.show()
    
    # Construct the instrument matrix
    A3 = np.array([(pd1R-pd1L)/2, (pd2R-pd2L)/2, (pd3R-pd3L)/2, (pd4R-pd4L)/2])
    A3_err=np.array([np.sqrt(pd1R_err**2+pd1L_err**2)/2, np.sqrt(pd2R_err**2+pd2L_err**2)/2,
                     np.sqrt(pd3R_err**2+pd3L_err**2)/2, np.sqrt(pd4R_err**2+pd4L_err**2)/2])
    
    A3 = np.matrix.transpose(np.array([A3]))
    A3_err = np.matrix.transpose(np.array([A3_err]))
    
    A02 = np.array(fit_parameters)  # data from linear calibration
    A02_err = np.array(fit_errs)
    
    A[:][:][q] = np.hstack((A02, A3))  # This is the instrument matrix!
    A_err = np.hstack((A02_err, A3_err))
    
 #   print('Instrument matrix A: ')
 #   print(A)
 #   print('')
    
#    A_cond=np.linalg.cond(A[:][:][q], p = 2)  # condition number with the standard 2-norm
#    print('Condition number of A: ')
#    print(A_cond)
#    print('')
    
    Ainv[:][:][q] = np.linalg.inv(A[:][:][q])  # this is the inverse of the instrument matrix
#    print('Inverted matrix Ainv: ')
#    print(Ainv[:][:][q])
#    print('')
    
    # save the instrument matrix as a text file for use in other scripts
    
    #if 'linux' in platform:
        #np.savetxt('../Ainv.txt', Ainv)
    #else:
        #np.savetxt('..\\Ainv.txt', Ainv)
    
#    np.savetxt('..\\Ainv.txt', Ainv[:][:][q])
    
    #Error in Ainv (see https://arxiv.org/pdf/hep-ex/9909031.pdf, http://sci-hub.io/10.1364/ao.47.002541)
    #Ainv_err=np.abs(np.dot(np.dot(Ainv, A_err),Ainv)) #need to change starting here
    
    #Assuming elements in A have no covariance between each other
#    Ainv_cov=np.zeros((4,4,4,4))
#    
#    # outer four sums to compute the whole covariance array
#    for aa in range(4):
#        for bb in range(4):
#            for a in range(4):
#                for b in range(4):
#                    # sum over i and j
#                    s=0.
#                    # inner two summations just for one element
#                    for i in range(4):
#                        for j in range(4):
#                            s += Ainv[aa][i]*Ainv[j][bb]*Ainv[a][i]*Ainv[j][b]*(A_err[i][j])**2
#                    Ainv_cov[aa][bb][a][b] = s
    
#    print('Covariance in Ainv: ')
#    print(Ainv_cov)
    #np.savetxt('..\\Ainv_cov.txt', Ainv_cov)
    # save the covariance matrix to a text-like file?
    #if 'linux' in platform:
    #    pickle.dump( Ainv_cov, open( "../Ainv_cov.p", "wb" ) )
    #else:
    #    pickle.dump( Ainv_cov, open( "..\Ainv_cov.p", "wb" ) )
    
#    pickle.dump( Ainv_cov, open( "..\Ainv_cov.p", "wb" ) )
    os.chdir('..')
    os.chdir('..')
    q+=1

B=np.zeros((16,len(angledirs)-1))    
for indeks in range(1,len(angledirs)):
 #   print(indeks)
    taeller=0
    for i1 in range(4):
        for i2 in range(4):
            B[taeller][indeks-1]=abs((A[indeks][i1][i2]-A[0][i1][i2])/A[0][i1][i2])
            taeller+=1
plt.figure()
z_min, z_max = 0,1# -np.abs(B).max(), np.abs(B).max()
plt.imshow(B,vmin=z_min, vmax=z_max)#,cmap='RdBu'
plt.colorbar()

print('Analyzer matrix deviation: ')
print(B)

#%%  plot
# plotting on Poincare sphere

#class Arrow3D(FancyArrowPatch):
#    def __init__(self, xs, ys, zs, *args, **kwargs):
#        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
#        self._verts3d = xs, ys, zs
#
#    def draw(self, renderer):
#        xs3d, ys3d, zs3d = self._verts3d
#        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
#        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
#        FancyArrowPatch.draw(self, renderer)
#
## Set the aspect ratio to 1 so our sphere looks spherical
#fig = plt.figure(figsize=plt.figaspect(1.))
#ax = fig.add_subplot(111, projection='3d')
#
#def plot_sphere(ax,arrows='xyz',equatorial=True):
#    phi = np.linspace(0, np.pi, 200)
#    theta = np.linspace(0, 2*np.pi, 200)
#
#    #equatorial circle
#    xe=np.sin(theta)
#    ye=np.cos(theta)
#
#    phi, theta = np.meshgrid(phi, theta)
#
#    # The Cartesian coordinates of the unit sphere
#    x = np.sin(phi) * np.cos(theta)
#    y = np.sin(phi) * np.sin(theta)
#    z = np.cos(phi)
#
#    ax.plot_surface(x, y, z,  rstride=10, cstride=10, color='#EBE3E8',
#                antialiased=True, alpha=0.5, lw=0.)#, facecolors=cm)
#    if 'y' in arrows:
#        ax.add_artist(Arrow3D([0, 0], [-0.03, 1.5], 
#                        [0,0], mutation_scale=15, 
#                        lw=0.25, arrowstyle="-|>", color="black"))
#        ax.text(0,1.5,0, '$S_2$', fontweight='bold')        
#    if 'x' in arrows:
#        ax.add_artist(Arrow3D([0.0, 1.5], [0,0], 
#                        [0,0], mutation_scale=15, 
#                        lw=0.25, arrowstyle="-|>", color="black"))
#        ax.text(1.6,0,0, '$S_1$', fontweight='bold')        
#    if 'z' in arrows:        
#        ax.add_artist(Arrow3D([0, 0], [0,0], 
#                        [-0.03,1.5], mutation_scale=15, 
#                        lw=0.25, arrowstyle="-|>", color="black"))
#        ax.text(0,0,1.5, '$S_3$',fontweight='bold')
#    if equatorial:
#        ax.plot(xe,ye,0,'--', dashes=(10, 10), lw=0.25, color='red', alpha=1)
#
#plot_sphere(ax)
#   
#for j in range(len(angledirs)):
#    for h in range(4):
#        ax.plot([A[1][h][j]/np.linalg.norm([A[1][h][j],A[][h][j],S3.item(j)])], [S2.item(j)/np.linalg.norm([S1.item(j),S2.item(j),S3.item(j)])], [S3.item(j)/np.linalg.norm([S1.item(j),S2.item(j),S3.item(j)])],  c = farve[i]+2*j/(3*len(S1)), marker='o'))#A[1:3][h][j])
#    ax.scatter(S1, S2, S3, marker='o')
##%% Define functions to reconstruct Stokes vector and compute DOP
#def determine_stokes(measurement):
#    try:
#        return np.dot(Ainv, measurement)
#    except ValueError:
#        raise('Input is not a 4 x 1 intensity measurement!')
#
#def determine_dop(pol_state):
#    return np.sqrt(pol_state[1]**2+pol_state[2]**2+pol_state[3]**2)/pol_state[0]
#        
##%% Now check to see that results of instrument matrix make sense by performing 
###  consistency check on linear states
#
#dops = []
#stokes_list = []
#for i in range(len(pd1_voltage)):
#    stokes = np.dot(Ainv, pd_voltages[:, i])
#    dop = np.sqrt(stokes[1]**2 + stokes[2]**2 + stokes[3]**2)/stokes[0]
#    dops.append(dop)
#    stokes_list.append(stokes/stokes[0])
#    
#plt.figure(2)
#plt.scatter(np.arange(len(dops)),dops)
#plt.plot([0,len(dops)],[1,1],alpha=0.3, color='black')
#plt.show()
#
#L = np.array([pd1L, pd2L, pd3L, pd4L])
#stokes=np.dot(Ainv, L)
#
