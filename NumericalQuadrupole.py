from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from nfft import nfft

#set up plotting environment
fig = plt.figure(figsize=(10,10))


#Define 3d/2d plotting


col = 'C0'
fs = 20
d = 3


ax1 = plt.subplot2grid((1,1), (0,0))



#PLOT MPD TRAJECTORY

MPDfile = 'TEST_INERTIA.txt'


data1 = np.loadtxt(MPDfile)
t = data1[:,0]
hp = data1[:,1]
hc = data1[:,2]

print (min(t),np.argmin(hp), min(hc))

ax1.plot(t,hp)
ax1.plot(t,hc)


plt.show()














