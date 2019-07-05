from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema
import scipy
from matplotlib.transforms import Bbox



plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#Set up plotting environment
fig = plt.figure(figsize=(20,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))





#Get Data
MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'

#Load data
data1 = np.loadtxt(MPDfile) #if we start at rp to set log of asc node, skip first bits
x = data1[:,0]
y = data1[:,1]
z = data1[:,2]

r = data1[:,12]
hPLUS = data1[:,23]
OBSR = data1[0,24]
convert_s = data1[0,25]


ind = np.argmin(r)
rmin = r[ind]
xmin = x[ind]
hmin = hPLUS[ind]

xran = np.arange(2, 1e7)
delta = scipy.integrate.simps(hmin/xran, xran)
print (delta/convert_s)

ax1.plot(xran, hmin/xran)
ax1.set_yscale('log')
plt.show()




