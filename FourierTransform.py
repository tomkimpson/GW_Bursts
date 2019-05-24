from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#SEt up plotting environment
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))

#Load data

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)


t = data1[:,3]
hplus = data1[:,4]
hcross = data1[:,5]


#ax1.plot(t,hplus, c='C0')
#ax1.plot(t,hcross, c='C1')


#Interpolate to get even sampling



Rate = 1

n = Rate*len(t)
t1 = np.linspace(t[0],t[-1],n)
hplus = np.interp(t1,t,hplus)
hcross= np.interp(t1,t,hcross)







dt = t1[1] - t1[0]

print ('dt = ', (dt))
print ('Sampling frequency = ', (1/dt))


#Calculate the FT
hplusT = dt*np.fft.rfft(hplus)
hcrossT = dt*np.fft.rfft(hcross)

#Get the frequencies
f = np.fft.rfftfreq(len(hplus), dt)
df = f[1] - f[0] #arethe frequencies evelyspaces?



#Get rid of zeroth frequencies

hplusT = hplusT[1:] # get rid of zeroth frequency
hcrossT = hcrossT[1:]
f = f[1:]

print ('Completed FT')

#ax1.scatter(t1, hplus,c='C0')
#ax1.scatter(t1, hcross,c='C1')


#LISA NOISE CURVE
Larm = 2.5e9
Clight = 3e8
fstar = Clight/(2*np.pi*Larm)
NC = 2

#Constants for glactic binary confusion noise
alpha = 0.138
beta = -221.
kappa = 521.
gamma = 1680.
f_knee = 1.13e-3

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




#Calculate the SNR


h = hplusT + hcrossT
hS = np.conj(h)


SNR2 = np.real(2*np.sum((hS*h + h*hS)/(S) * df))

SNR = np.sqrt(SNR2)

print ('The calculated SNR = ', SNR)



'''
ax2.scatter(f,htilde1)
ax2.scatter(f,htilde2)
ax2.loglog(f,np.sqrt(f*S))
ax2.set_xlim(1.0e-6,1.0e0)
ax2.set_ylim(1e-22,1.0e-13)
'''









