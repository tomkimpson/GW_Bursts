from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema
import scipy
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#SEt up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))
#ax2 = plt.subplot2grid((1,2), (0,1))


#Load data
MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'
data1 = np.loadtxt(MPDfile)
t = data1[:,3]
hplus = data1[:,4]
hcross = data1[:,5]
r = data1[:,12]
phi = data1[:,14]
factor = data1[0,15]
#Get the analytical model 
hpA = -4/r[0] *np.cos(2*phi)
hcA = -4/r[0] *np.sin(2*phi)
x = data1[:,0]
y = data1[:,1]
z = data1[:,2]



'''
#Plot the trajectory
ax1.plot(x,y)
ax1.scatter(0,0,c='r')

'''

#Plot the waveform
hpNORMp = data1[:,6]
hpNORMc = data1[:,7]


ax1.plot(t/(24*3600),hpNORMp, c='C0')
#ax2.plot(t/(24*3600),hpNORMc, c='C1')
ax1.plot(t/(24*3600), hpA, c='C2')
#ax2.plot(t/(24*3600), hcA, c='C3')






#Now calculate the overlap

Tobs = 1*24*60*60 # 1 day observation
#and scale properly

hpA = hpA /factor 
hcA = hcA /factor 






#Interpolate to get even sampling
fs = 2 #sampling frequency
dt = 1.0/fs
maxF = 1.0
t1 = np.arange(t[0], Tobs, dt)
hplus = np.interp(t1,t,hplus)
hcross= np.interp(t1,t,hcross)
hpA = np.interp(t1,t,hpA)
hcA = np.interp(t1,t,hcA)





#Get the frequencies
f = np.fft.rfftfreq(hplus.size, dt)
df = f[1] - f[0] #arethe frequencies evelyspaces?


#Calculate the FT
hplusT = dt*np.fft.rfft(hplus) #/ factorW
hcrossT = dt*np.fft.rfft(hcross) #/ factorW

hpA_T = dt*np.fft.rfft(hpA) #/ factorW
hcA_T = dt*np.fft.rfft(hcA) #/ factorW

#Get rid of zeroth frequencies

hplusT = hplusT[1:] # get rid of zeroth frequency
hcrossT = hcrossT[1:]
hpA_T = hpA_T[1:]
hcA_T = hcA_T[1:]



f = f[1:]
print ('Completed FT')





#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


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
RFILE = np.loadtxt('/unsafe/tok2/ResponseFunction.txt')
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


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------





hNumerical = np.sqrt(abs(hplusT)**2 + abs(hcrossT)**2)
hAnalytical = np.sqrt(abs(hpA_T)**2 + abs(hcA_T)**2)




#overlap = 2 * scipy.integrate.simps( (hNumerical * np.conj(hAnalytical) + np.conj(hNumerical) * hAnalytical)/(S),f)


normN = 2 * scipy.integrate.simps( (hNumerical * np.conj(hNumerical) + np.conj(hNumerical) * hNumerical)/(S),f)
normA = 2 * scipy.integrate.simps( (hNumerical * np.conj(hNumerical) + np.conj(hNumerical) * hNumerical)/(S),f)



hNumerical = hNumerical/np.sqrt(normN)
hAnalytical = hAnalytical/np.sqrt(normN)



overlap = 2 * scipy.integrate.simps( (hNumerical * np.conj(hAnalytical) + np.conj(hNumerical) * hAnalytical)/(S),f)


print ('overlap:', overlap)

#Make it look pretty

fs=20

'''
ax1.set_xlabel(r'$ x [\rm r_g]$',fontsize = fs)
ax1.set_ylabel(r'$ y [\rm r_g]$',fontsize = fs)
ax1.locator_params(axis='both', nbins = 8)
ax1.tick_params(axis='both', which='major', labelsize=16)
'''

ax1.set_xlim(-0.002,0.1)
ax1.set_xlabel(r' t [days]',fontsize = fs)
ax1.set_ylabel(r'$ h_{+} (r/\mu)$',fontsize = fs)
ax1.locator_params(axis='both', nbins = 8)
ax1.tick_params(axis='both', which='major', labelsize=16)


plt.savefig('/unsafe/tok2/PAPERS/GW_Burst/overlap_waveform.png', dpi=300)



plt.show()




