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


#Set up plotting environment
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))
#Load data

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


#Load the data and plot the waveform
data1 = np.loadtxt(MPDfile)
t = data1[:,3]

#theta = 0
hplus1 = data1[:,6]
hcross1 = data1[:,7]


#theta = pi/4
hplus2 = data1[:,8]
hcross2 = data1[:,9]



#theta = pi/2
hplus3 = data1[:,10]
hcross3 = data1[:,11]




ax1.plot(t,hplus1, c='C0')
ax1.plot(t,hcross1, c='C1')




plt.show()




