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

#Set the path
path = os.environ['GWDir']
files=glob.glob(path+'*.txt')


def waveform(f):
    data = np.loadtxt(f)
    
    t = data[:,0]
    hp = data[:,4]
    hc = data[:,5]

    ax.plot(t,hp)
    ax.plot(t,hc)




for f in files:
   waveform(f)


plt.show()


