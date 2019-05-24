from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema



#SEt up plotting environment
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))

#Load data

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)

hplus = data1[:,3] # this is for obs theta = 0
hcross = data1[:,4]
t = data1[:,9]


ax1.plot(t,hplus)





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


ax2.scatter(freq,htilde)



plt.show()




