from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os





d = 2


#Set up plotting environment
if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0))



#Load data


path = os.environ["GWDir"]



fileA = path+'ForPaper/Fig3a.txt'

data = np.loadtxt(fileA)

t = data[:,0]
x = data[:,1]
y = data[:,2]
z = data[:,3]



#Plot it


if (d == 3):
    ax1.plot(x,y,z)  
    ax1.scatter(x[0],y[0],z[0],c='r')  
    limit = max(max(x),max(y),max(z))
    ax1.set_xlim(-limit,+limit)
    ax1.set_ylim(-limit,+limit)
    ax1.set_zlim(-limit,+limit)

if (d == 2):
    ax1.plot(x,y)
    ax1.scatter(0,0,c='r')
    lower = -1600
    upper = 500
    ax1.set_xlim(lower,upper)
    ax1.set_ylim(lower,upper)





plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20
ax1.set_xlabel(r'$x [r_g]$', fontsize = fs)
ax1.set_ylabel(r'$y [r_g]$', fontsize = fs)
ax1.locator_params(axis='both', nbins=5)
ax1.tick_params(axis='both', which='major', labelsize=fs-4)

plt.savefig('/Users/tomkimpson/Dropbox/MSSL/Papers/PaperNGW_burst/figures/Fig3Trajectory.png',dpi=300)
plt.show()
