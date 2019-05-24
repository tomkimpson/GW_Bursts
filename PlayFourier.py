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













N = 8
dt = 0.1
df = 1/(dt*N)
a = np.random.normal(size=N)
b = np.fft.fft(a)

E1 = np.sum(a*np.conj(a)*dt)
E2 = np.sum(b*np.conj(b)*df)


print (E1, E2)





