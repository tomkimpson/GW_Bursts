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
ax1 = plt.subplot2grid((3,1), (0,0))
ax2 = plt.subplot2grid((3,1), (1,0),sharex=ax1)
ax3 = plt.subplot2grid((3,1), (2,0),sharex=ax1)



#PLOT MPD TRAJECTORY

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)



t = data1[:,3] / (60*60*24*365)

#theta = 0
hplus1 = data1[:,6]
hcross1 = data1[:,7]


#theta = pi/4
hplus2 = data1[:,8]
hcross2 = data1[:,9]



#theta = pi/2
hplus3 = data1[:,10]
hcross3 = data1[:,11]



#plot it all
ax1.plot(t,hplus1, c='C0')
ax1.plot(t,hcross1, c='C1')


ax2.plot(t,hplus2, c='C0')
ax2.plot(t,hcross2, c='C1')


ax3.plot(t,hplus3, c='C0')
ax3.plot(t,hcross3, c='C1')





#Sort teh axes

ax3.set_xlabel(r' t [years]',fontsize = fs)


ax1.set_ylabel(r'$ h_{+, \times} (r/\mu)$',fontsize = fs)
ax2.set_ylabel(r'$ h_{+, \times} (r/\mu)$',fontsize = fs)
ax3.set_ylabel(r'$ h_{+, \times} (r/\mu)$',fontsize = fs)


plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

ax1.locator_params(axis='y', nbins = 3)
ax2.locator_params(axis='y', nbins = 3)
ax3.locator_params(axis='y', nbins = 3)

ax3.locator_params(axis='x', nbins = 5)

ax1.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax3.tick_params(axis='both', which='major', labelsize=16)


limit = 1.0
ax1.set_ylim(-limit,limit)
ax2.set_ylim(-limit,limit)
ax3.set_ylim(-limit,limit)



plt.savefig('/unsafe/tok2/PAPERS/GW_Burst/Waveform.png', dpi=300)




plt.show()




