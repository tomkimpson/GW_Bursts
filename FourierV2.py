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


#Set up plotting environment
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))

#Load data

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)


t = data1[:,3]

hplus = np.cos(t)
#hplus = data1[:,4]
hcross = np.sin(t)






t0 = -100
dt = 0.001
t = np.arange(t0,-t0,dt)

f = 1/(t**2 + 1)

ax1.plot(t,f)




#Now get FT
#g = np.fft.fft(f)
#w = np.fft.fftfreq(len(f)*2*np.pi/dt)





plt.show()













def Fourier(t,h):

    #Interpolate to get even sampling
    n = 10*len(t)
    newt = np.linspace(t[0],t[-1],n)
    newh = np.interp(newt,t,hplus)

    t1 = newt
    h = newh 
    dt = t1[1] - t1[0]

   #Calculate the Fourier transform
    htilde = dt*np.fft.rfft(h)
    freq = np.fft.rfftfreq(len(h),dt)

    htilde = htilde[1:] # get rid of zeroth frequency
    freq = freq[1:]
    return t1,h,freq,htilde


t1,h1,freq1,htilde1 = Fourier(t,hplus)
t2,h2,freq2,htilde2 = Fourier(t,hcross)



ax1.plot(freq1,htilde1)
plt.show()




sys.exit()



H = htilde1 + htilde2 # h(f)
HStar = np.conj(H)
hmag = H*HStar


f = freq1
logf = np.log(f)



df = f[1] - f[0] #approx
d_logf = logf[1] - logf[0]




#Plot the waveform in teh time and frequency domains

ax1.plot(t1,h1)
ax1.plot(t2,h2)



#Calculate the SNR
hstar1 = np.conj(htilde1)
hstar2 = np.conj(htilde2)



ax2.axvline(freq1[-1])

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



#Barackandcutler noise
'''
Sinst = 1.22e-51 * f**(-4) + 2.12e-41 + 1.22e-37 * f**2
Sgal = 2.1e-45 *(f)**(-7/3)

SEG =4.2e-47 * (f)**(-7/3)


dNdf = 2e-3 *(1/f)**11/3

print (f[0])

print (dNdf[0])

print (f[0]**(-11/3) )

sys.exit()

print (np.exp(-1.5*dNdf))

SIG = np.empty_like(f)
for i in range(len(f)):
    SIG[i]= min(Sinst[i]/np.exp(-1.5*dNdf[i]),Sinst[i] + Sgal[i])



S = SIG + SEG
'''








ax2.scatter(f,htilde1)
ax2.scatter(f,htilde2)
ax2.loglog(f,np.sqrt(f*S))
ax2.set_xlim(1.0e-6,1.0e0)
ax2.set_ylim(1e-22,1.0e-13)



#SNR2 = 4*np.sum(f*d_logf*hmag/Pn)


SNR2 =4*np.real(np.trapz(df*hmag/Pn))




SNR = np.sqrt(SNR2)

print (SNR)






#hmag = abs(htilde1)**2
#SNR2 = 4*np.sum(hmag/Pn)
#print (SNR2, np.sqrt(SNR2))


#SNR2 = 4*np.trapz(hmag,abs(Pn))
#print (np.sqrt(SNR2))




#hmag = htilde1*hstar1 + htilde2*hstar2


#SNR2 = 4*np.real(np.trapz(hmag/S))

#print ('SNR', np.sqrt(SNR2))


plt.show()









