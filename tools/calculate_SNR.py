from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sys


#Load the data
path = os.environ["GWDir"]
files = glob.glob(path+'*.txt')


#Set up plotting environment
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))





def process(f):
    
    #Load the data
    data = np.loadtxt(f)
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
    Tobs = 1*24*60*60 # 1 day observation
    ax1.axvline((tmin+Tobs/2), linestyle = '--', c='0.5')
    ax1.axvline((tmin-Tobs/2), linestyle = '--', c='0.5')


    #The data was generated using an adaptive stepsize method
    #Therefore interpolate to get even sampling 

    fs = 2.0 #sampling frequency
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
    print ('LENGTHS:', len(hplus))

    #Calculate the FT
    hplusT = dt*np.fft.rfft(hplus) #/ factorW
    hcrossT = dt*np.fft.rfft(hcross) #/ factorW
    
    #The net signal
    hsig = abs(hplusT)**2 + abs(hcrossT)**2
    ax2.loglog(f,np.sqrt(hsig), C = 'C3')


    print (hcrossT)
    print ('Got the FT')
for f in files:
    process(f)




plt.show()    
