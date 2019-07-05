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
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))


#Load data
MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'
data1 = np.loadtxt(MPDfile)
t = data1[:,3]
hplus = data1[:,4]
hcross = data1[:,5]
r = data1[:,12]

#Plot the waveform

norms = hplus[0]

hpNORM = data1[:,6]
ax1.plot(t,hplus/hplus[0], c='C0')
ax1.axhline(4/r[0], linestyle = '--',c='0.5')


#Only get the data within Tobs/2 of periapsis

ind = np.argmin(r)
rmin = r[ind]
tmin = t[ind]
Tobs = 1*24*60*60 # 1 day observation



ax1.axvline((tmin+Tobs/2), linestyle = '--', c='0.5')
ax1.axvline((tmin-Tobs/2), linestyle = '--', c='0.5')








#Interpolate to get even sampling
fs = 2.0 #sampling frequency
dt = 1.0/fs
t1 = np.arange(tmin-Tobs/2, tmin+Tobs/2, dt)
hplus = np.interp(t1,t,hplus)
hcross= np.interp(t1,t,hcross)

ax1.scatter(t1, hplus/norms,c='C1')



#Get the frequencies
f = np.fft.rfftfreq(hplus.size, dt)
df = f[1] - f[0] #arethe frequencies evelyspaces?
print ('LENGTHS:', len(hplus))






#Calculate the FT
hplusT = dt*np.fft.rfft(hplus) #/ factorW
hcrossT = dt*np.fft.rfft(hcross) #/ factorW






print ('Got the FT')






#Get rid of zeroth frequencies

hplusT = hplusT[1:] # get rid of zeroth frequency
hcrossT = hcrossT[1:]
f = f[1:]
print ('Completed FT')



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




#R = 0.30/(1+0.6*(f/fstar)**2)











#Power Spectral density
P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2



#Total noise
S = Pn/R + Sc

#S = 10/(3*Larm**2) * (P_oms + 2*(1+np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4) * (1+0.60*(f/fstar)**2)






SNR2 = 4 * scipy.integrate.simps((abs(hplusT)**2 + abs(hcrossT)**2)/S , f)
SNR = np.sqrt(SNR2)
print ('The calculated SNR = ', SNR)


#ax2.loglog(f, Pn/R)

#ax2.scatter(f, (abs(hplusT)**2 + abs(hcrossT)**2)/S)

#plt.show()

ax2.loglog(f,abs(hplusT))
ax2.scatter(f,abs(hplusT))


g = open('SNR_Data.txt', 'a')
g.write(str(r[0]))
g.write(' ')
g.write(str(SNR))
g.write('\n ')
g.close()


ax2.loglog(f,np.sqrt(f*S))
ax2.loglog(f,np.sqrt(S))
ax2.set_xlim(1.0e-6,1.0e0)
ax2.set_ylim(1e-22,1.0e-13)




plt.show()




