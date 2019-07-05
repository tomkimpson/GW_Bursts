from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from nfft import nfft






plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

#Set up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))



f = np.logspace(np.log10(1.0e-5), np.log10(1.0e0),1000)



#LISA NOISE CURVE
Larm = 2.5e9
Clight = 3e8
fstar = Clight/(2*np.pi*Larm)
NC = 2

#Constants for glactic binary confusion noise
alpha = 0.133
beta = 243.
kappa = 482.
gamma = 917.
f_knee = 2.58e-3

A = 1.8e-44/NC
Sc = 1. + np.tanh(gamma*(f_knee-f))
Sc *=np.exp(-f**alpha + beta*f*np.sin(kappa*f))

Sc *= A*f**(-7./3.)


#Response function load from folder https://github.com/eXtremeGravityInstitute/LISA_Sensitivity 
RFILE = np.loadtxt('ResponseFunction.txt')
Rx = RFILE[:,0] * fstar
Ry = RFILE[:,1] * NC

newR = np.interp(f,Rx,Ry)
R = newR



#Power Spectral density
P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2



#Total noise
S = Pn/R + Sc

#S = 10/(3*Larm**2) * (P_oms + 2*(1+np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4) * (1+0.60*(f/fstar)**2)


ax1.loglog(f,np.sqrt(f*S))
ax1.set_xlim(1.0e-5,1.0e0)
ax1.set_ylim(3.0e-22,1.0e-15)







ax1.set_xlabel(r' f [Hz]',fontsize = fs)
ax1.set_ylabel(r'Characteristic Strain',fontsize = fs)
ax1.tick_params(axis='both', which='major', labelsize=16)

plt.savefig('/unsafe/tok0/PAPERS/GW_Burst/LISANoise.png', dpi=300)

plt.show()






