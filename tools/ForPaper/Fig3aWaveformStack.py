from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
import os

#Set up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((3,1), (0,0))
ax2 = plt.subplot2grid((3,1), (1,0),sharex=ax1)
ax3 = plt.subplot2grid((3,1), (2,0),sharex=ax1)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20



#Set the path
path = os.environ['GWDir']
Fig = 'Fig2b'


datafile1=path+'ForPaper/Fig3a.txt'
datafile2=path+'ForPaper/Fig3b.txt'
datafile3=path+'ForPaper/Fig3c.txt'


def waveform(f,ax):
    
    #Load the data
    #length = np.loadtxt(f)
    #length = len(length[:,0])

    #data = np.loadtxt(f, skiprows=int(length/4)) #skip the start since we start at rp
    data = np.loadtxt(f)
    t = data[:,0]
    hp = data[:,4]
    hc = data[:,5]
    
    
    #Normalise to t=0
    r = data[:,8]
    ind = np.argmin(r)
    rmin = r[ind]
    tmin = t[ind]
    t = t-tmin

    #Convert from seconds to days
    t = t/(60*60*24)

    #Plot it
    ax.plot(t,hp)
    ax.plot(t,hc)

    #Make it pretty
    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both', which='major', labelsize=fs-4)
    ax.set_ylabel(r'$h_{+, \times} (r/\mu)$', fontsize = fs)
    Tobs = 50 #day
    #ax.set_xlim(-Tobs, +Tobs)
    ax.set_ylim(-0.070, +0.05)




waveform(datafile1,ax1)
waveform(datafile2,ax2)
waveform(datafile3,ax3)


plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax2.get_xticklabels(),visible=False)
ax3.set_xlabel('t [days]', fontsize = fs)
plt.savefig('/Users/tomkimpson/Dropbox/MSSL/Papers/PaperNGW_burst/figures/Fig3_stack.png',dpi=300)
plt.show()


