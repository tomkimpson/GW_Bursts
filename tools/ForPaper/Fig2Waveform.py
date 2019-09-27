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
ax = plt.subplot2grid((1,1), (0,0))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')




#Set the path
path = os.environ['GWDir']
Fig = 'Fig2b'


datafile=path+'ForPaper/'+Fig+'.txt'


def waveform(f):
    
    #Load the data
    length = np.loadtxt(f)
    length = len(length[:,0])

    data = np.loadtxt(f, skiprows=int(length/4)) #skip the start since we start at rp
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
    fs = 20
    ax.set_xlabel('t [days]', fontsize = fs)
    ax.set_ylabel(r'$h_{+, \times} (r/\mu)$', fontsize = fs)
    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both', which='major', labelsize=fs-4)

    Tobs = 5 #day
    ax.set_xlim(-Tobs, +Tobs)
    ax.set_ylim(-0.070, +0.05)




waveform(datafile)

plt.savefig('/Users/tomkimpson/Dropbox/MSSL/Papers/PaperNGW_burst/figures/'+Fig+'.png',dpi=300)
plt.show()


