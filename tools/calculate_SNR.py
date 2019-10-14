from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sys
import scipy.integrate

#Load the data
path = os.environ["GWDir"]
files = glob.glob(path+'*.txt')


#Set up plotting environment
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))





def process(f,Tintegration):
 
    #ignore the first burst - useful if we start at periapsis
    length = np.loadtxt(f)
    length = len(length[:,0])


    #Load the data
    data = np.loadtxt(f,skiprows=int(length/4))
    #data = np.loadtxt(f)
    t = data[:,0]
    hplus_norm = data[:,4]
    hcross_norm = data[:,5]
    hplus = data[:,6]
    hcross = data[:,7]
    r = data[:,8]
    N = data[0,9]




    #plot the waveform
    ax1.plot(t, hplus_norm)
    ax1.plot(t, hcross_norm)

    #Only select some of the data
    #that within a certain limit of periapsis

    ind = np.argmin(r)
    rmin = r[ind]
    tmin = t[ind]
    Tobs = Tintegration*24*60*60 # 1 day observation
    ax1.axvline((tmin+Tobs/2), linestyle = '--', c='0.5')
    ax1.axvline((tmin-Tobs/2), linestyle = '--', c='0.5')


    if (tmin+Tobs/2 > t[-1]):
        print ('Your time series does not extend that far. Integrate for longer')
        sys.exit()



    #The data was generated using an adaptive stepsize method
    #Therefore interpolate to get even sampling 

    fs = 2**1 #sampling frequency
    dt = 1/fs
    t1 = np.arange(tmin-Tobs/2, tmin+Tobs/2, dt)
    hplus = np.interp(t1,t,hplus)
    hcross= np.interp(t1,t,hcross)

    #if you want to check the interpolation graphiclly:
    #ax1.plot(t1,hplus/N)
    #ax1.plot(t1,hcross/N)


    #Get the frequencies
    f = np.fft.rfftfreq(hplus.size, dt)
    df = f[1] - f[0] #arethe frequencies evelyspaces?
    

    #Calculate the FT
    hplusT = dt*np.fft.rfft(hplus) #/ factorW
    hcrossT = dt*np.fft.rfft(hcross) #/ factorW
 
    
    #Get rid of zeroth frequencies - WHY?
    hplusT = hplusT[1:] # get rid of zeroth frequency
    hcrossT = hcrossT[1:]
    f = f[1:]


    #Calculate the LISA noise curve
    Larm = 2.5e9
    Clight = 3e8
    fstar = Clight/(2*np.pi*Larm)
    NC = 2



    #LISA response function
    RFILE = np.loadtxt('noise/ResponseFunction.txt')
    Rx = RFILE[:,0] * fstar
    Ry = RFILE[:,1] * NC


    #Get rid of all frequencies greater than Rx[-1]

    hplusT_temp = []
    hcrossT_temp = []
    f_temp = []
    for i in range(len(f)):
        if f[i] <= Rx[-1]:
            hplusT_temp.extend([hplusT[i]])    
            hcrossT_temp.extend([hcrossT[i]])    
            f_temp.extend([f[i]])    
        
    hplusT = np.array(hplusT_temp)
    hcrossT = np.array(hcrossT_temp)
    f = np.array(f_temp)










    alpha = 0.133
    beta = 243.
    kappa = 482.
    gamma = 917.
    f_knee = 2.58e-3

    A = 1.8e-44/NC
    Sc = 1. + np.tanh(gamma*(f_knee-f))
    Sc *=np.exp(-f**alpha + beta*f*np.sin(kappa*f))

    Sc *= A*f**(-7./3.)




    #Get rid of frequencies higher than limit set by Rx[-1]



 #   print (Rx[0], Rx[-1])
  #  sys.exit()

    newR = np.interp(f,Rx,Ry)
    R = newR



    #Power Spectral Density
    P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
    P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
    Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2

    #Total noise
    S = Pn/R + Sc



    #The net signal
    hsig = abs(hplusT)**2 + abs(hcrossT)**2

    SNR2 = 4 * scipy.integrate.simps((hsig)/S , f)
    SNR = np.sqrt(SNR2)
    print ('The calculated SNR = ', SNR)



    #some plotting
    ax2.loglog(f,np.sqrt(hsig), C = 'C3')
    ax2.loglog(f,np.sqrt(S), C = 'C3')


 #   print (hcrossT)
  #  print ('Got the FT')
Tint = 0.001
print (1/Tint)
for f in files:
    process(f, Tint)




plt.show()    
