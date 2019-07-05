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



#Set observer location 0, pi/4, pi/2
BigTheta = np.pi/4 #observer latitude


#Get Data
MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'
length = np.loadtxt(MPDfile)
length = len(length[:,0])

#Load data
data1 = np.loadtxt(MPDfile, skiprows = int(length/2)) #if we start at rp to set log of asc node, skip first bits

t = data1[:,3]


#Normalise to t=0 and only get the data within Tobs/2 of periapsis

r = data1[:,12]
ind = np.argmin(r)
Eloss = data1[:,17]
EK = data1[:,18]
Lx = data1[:,19]
Ly = data1[:,20]
Lz = data1[:,21]
LK = data1[:,22]




Tobs = 1.00*24*60*60 

#Get the data witin the bounds to use for SNR calculations

tmin = t[ind]
t_upper = tmin + Tobs/2
t_lower = tmin - Tobs/2

tt = []
EDOT = []
LxC = []
LyC = []
LzC = []
for i in range(len(t)):
    if t_lower<t[i]<t_upper:
        tt.extend([t[i]])
        EDOT.extend([Eloss[i]])
        LxC.extend([Lx[i]])
        LyC.extend([Ly[i]])
        LzC.extend([Lz[i]])





avX = np.average(LxC)
avY = np.average(LyC)
avZ = np.average(LzC)

Ldot = np.sqrt(avX**2 + avY**2 + avZ**2)

#print (np.average(LxC), np.average(LyC), np.average(LzC))


print (np.average(EDOT),EK[0])
print (avZ, LK[0])




ax1.plot(tt,EDOT)
ax2.plot(tt,LzC)
plt.show()


