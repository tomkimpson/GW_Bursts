from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sys
import scipy.integrate

#Load the data
path = os.environ["GWDir"]
datafile = path+'ForPaper/circular.txt'


#Set up plotting environment
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))




def FT(t, t1, h,dt):
    hnew = np.interp(t1,t,h)

    hT = dt*np.fft.rfft(hnew) #Transposed

    hT = hT[1:]


    return hT




def process(f):
    
    #Load the data
    data = np.loadtxt(f)
    t = data[:,0]
    hplus_norm = data[:,4]
    hcross_norm = data[:,5]
    r = data[:,8]
    N = data[0,9]
    phi = data[:,10]

    hplus_analytical = - 4/r[0] * np.cos(2*phi)
    hcross_analytical = - 4/r[0] * np.sin(2*phi)

    fs = 2.0
    dt = 1.0/fs
    t1 = np.arange(t[0], t[-1], dt)

       
    #Convert from seconds to hours
    tplot = t/(60*60)
    ax1.plot(tplot,hplus_analytical)    
    ax1.plot(tplot,hplus_norm)    
    
    f = np.fft.rfftfreq(t1.size,dt)
    f = f[1:]




    #Fourier transform the numerical and analytical signals
    hplusN = FT(t,t1,hplus_norm,dt)
    hcrossN = FT(t,t1,hcross_norm,dt)
    hplusA = FT(t,t1,hplus_analytical,dt)
    hcrossA = FT(t,t1,hcross_analytical,dt)


    hN = np.sqrt(abs(hplusN)**2 + abs(hcrossN)**2) #numerical
    hA = np.sqrt(abs(hplusA)**2 + abs(hcrossA)**2) #analytical








    #Noise

    #Calculate the LISA noise curve
    Larm = 2.5e9
    Clight = 3e8
    fstar = Clight/(2*np.pi*Larm)
    NC = 2

    alpha = 0.133
    beta = 243.
    kappa = 482.
    gamma = 917.
    f_knee = 2.58e-3

    A = 1.8e-44/NC
    Sc = 1. + np.tanh(gamma*(f_knee-f))
    Sc *=np.exp(-f**alpha + beta*f*np.sin(kappa*f))

    Sc *= A*f**(-7./3.)


    #LISA response function

    RFILE = np.loadtxt('../noise/ResponseFunction.txt')
    Rx = RFILE[:,0] * fstar
    Ry = RFILE[:,1] * NC

    newR = np.interp(f,Rx,Ry)
    R = newR


    #Power Spectral Density
    P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
    P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
    Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2

    #Total noise
    S = Pn/R + Sc
    S = Pn


    normN = 2 * scipy.integrate.simps( (hN * np.conj(hN) + np.conj(hN) * hN)/(S),f)
    normA = 2 * scipy.integrate.simps( (hA * np.conj(hA) + np.conj(hA) * hA)/(S),f)





    hN = hN / np.sqrt(normN)
    hA = hA / np.sqrt(normA)

    
    overlap = 2 * scipy.integrate.simps( (hN * np.conj(hA) + np.conj(hN) * hA)/(S),f)


    print ('The overlap is ', overlap)




process(datafile)

fs = 20
ax1.set_xlabel('t [hours]', fontsize = fs)
ax1.set_ylabel(r'$h_{+} (r/\mu)$', fontsize = fs)
ax1.locator_params(axis='both', nbins=5)
ax1.tick_params(axis='both', which='major', labelsize=fs-4)

plt.savefig('/Users/tomkimpson/Dropbox/MSSL/Papers/PaperNGW_burst/figures/overlap.png',dpi=300)
plt.show()    
