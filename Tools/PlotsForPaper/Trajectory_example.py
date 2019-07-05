from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from nfft import nfft






plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

#Set up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))



#PLOT MPD TRAJECTORY

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)
x = data1[:,0]
y = data1[:,1]
z = data1[:,2]


ax1.plot(x,y,c='C0')
#ax1.scatter(x[0],y[0],c='C4')
    
    
ax1.set_xlabel(r'$ x [\rm r_g]$',fontsize = fs)
ax1.set_ylabel(r'$ y [\rm r_g]$',fontsize = fs)




ax1.locator_params(axis='both', nbins = 5)

ax1.tick_params(axis='both', which='major', labelsize=16)



ax1.scatter(0,0,c='r')




plt.savefig('/unsafe/tok2/PAPERS/GW_Burst/Trajectory.png', dpi=300)


plt.show()








plt.show()




