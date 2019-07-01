# -*- coding: utf-8 -*-
"""
Spyder Editor

"""

import matplotlib.pyplot as plt
import numpy as np
import os
import tkinter as tk
from scipy.optimize import curve_fit
from scipy.special import erf
import math
from scipy.signal import chirp, find_peaks, peak_widths

plt.close()
print("*******SPOT SIZER*******\nBy Mauricio Bejarano\n",end='\n')

#get user inputs, directory, file output name
go_to_dir=input("Please provide directory path where files are located: ")
#go to directory containing files
#os.chdir("/Users/mbejarano/Documents/M2_express_linescanx_sweepz")

x_name=input("Please specify name of file with x values: ")
y_name=input("Please specify name of file with y values: ")
os.chdir(go_to_dir)
output_file_name=input("Please provide name for output file:(no extension) ")
coordinate=input("Please indicate whether it was a scan in x or in y: ")

#Extract X coordinates
X = np.genfromtxt(f'{go_to_dir}/{x_name}',usecols=1)
#X=raw_data_dimx[:,1]

#Read y values
Y=np.genfromtxt(f'{go_to_dir}/{y_name}',usecols=1)
#Y=raw_data_dimy[:,1]
#plt.plot(X, Y,'o',label='data')

#Functions for fitting curve
def origin_funct(x, A1,A2,L,p):
    return A1 + (A2-A1)/(1 + 10**((L-x)*p))
def gauss(x,y0,xc,w,A):
    return y0 + (A/(w*math.sqrt(math.pi/2)))*math.exp(-2*((x-xc)/w)**2)
gauss2 = np.vectorize(gauss)

#Scaling function
def scale(X, x_min, x_max):
    nom = (X-X.min())*(x_max-x_min)
    denom = X.max() - X.min()
    denom = denom +(denom is 0)
    return x_min + nom/denom 

# Fit the curve
params2, extras2 = curve_fit(origin_funct, X, Y, method='lm',p0=[0.09,0.24,5560,1.34])
#p0=[0.2,0.3,71919,2.68]

# Calculate the fitting
y_fitted4 = origin_funct(X, *params2)


#Calculate the derivative of the fitting curve
deriv=np.gradient(y_fitted4)


#Actual plotting
fig=plt.figure()
plt.clf()
ax1=plt.gca()
plt.xlabel(f'{coordinate} position ($\mu$m)')
plt.ylabel('Intensity (arb.u.)')
ax1.title.set_text('Knife edge scans')
ax1.plot(X, scale(Y,0,1),'o',label='Measured data',color='k') #Raw data
ax1.plot(X,scale(y_fitted4,0,1),'-', label='Sigmoidal fit',color='k', linewidth=0.7)      #Sigmoidal fitting


#ax2.plot(X,scale(deriv,-1,1),'-', linewidth=4)           #Derivative of the fitting
ax1.plot(X,scale(deriv,0,1),'-', linewidth=1, label='Derivative of Gaussian error function', color='b')
#scale(deriv,0,1)
#Fit derivative to Gaussian function
#params3, extras3 = curve_fit(gauss2, X, deriv, method='lm',p=[9.731E-4,5560,1.015,0.144])
#p0=[2.163E-4,71919,0.6,0.09]

#Calculate fitted gaussian curve
#gauss_fitted = gauss2(X,*params3)

#Plot fitted Gaussian bell
#ax2.plot(X,gauss_fitted,'-.',linewidth=2)

#Calculate width in samples of FWHM
peak, _ = find_peaks(deriv)
results_half = peak_widths(deriv, peak, rel_height=0.5)
#print(results_half[0])  # widths

#Calculate FWHM
delta_x=X[1]-X[0]
fwhm=results_half[0]*delta_x
limit_left=X[int(results_half[2][0])]
limit_right=X[int(results_half[3][0])]
#plt.hlines(0.5, limit_left, limit_right, color="b",)
print(f'FWHM is: {fwhm[0]:.3f} \u03BCm\n')
print('Graph saved at: ',go_to_dir,end='\n')
ax1.text(0.1, 0.6, f'FWHM={fwhm[0]:.3f} $\mu$m', transform=ax1.transAxes)
ax1.legend()
#fig.tight_layout()
plt.savefig(f'{go_to_dir}/{output_file_name}.png',bbox_inches = 'tight')
#fig.show()

