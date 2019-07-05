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


#Set observer location 0, pi/4, pi/2
BigTheta = 0 #observer latitude

#Get Data
MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'
length = np.loadtxt(MPDfile)
length = len(length[:,0])

data1 = np.loadtxt(MPDfile, skiprows = int(length/2)) #if we start at rp to set log of asc node, skip first bits
t = data1[:,3] / (60*60*24) #days

x = data1[:,0]
y = data1[:,1]
z = data1[:,2]
r = data1[:,12]



if BigTheta == 0:
    hplus1 = data1[:,6]
    hcross1 = data1[:,7]
elif BigTheta == np.pi/4:
    hplus1 = data1[:,8]
    hcross1 = data1[:,9]
elif BigTheta == np.pi/2:
    hplus1 = data1[:,10]
    hcross1 = data1[:,11]


#Normalise to t=0
index = np.argmin(r)
t0 = t[index]
t = t- t0





#plot it all
ax1.plot(t,hplus1, c='C0')
ax1.plot(t,hcross1, c='C1')


#Calculate the inclination w.r.t observer
a = np.array([x[0], y[0], z[0]])
b = np.array([x[45], y[45], z[45]])
Lvector = np.cross(a,b)
Lvector = Lvector/np.linalg.norm(Lvector)


Nvector = np.array([np.sin(BigTheta), np.sin(BigTheta), np.cos(BigTheta)])

inc = np.arccos(np.dot(Lvector, Nvector))#
print ('Inclination w.r.t observer = ', inc * 180/np.pi)







ax1.set_xlabel(r' t [days]',fontsize = fs)
ax1.set_ylabel(r'$ h_{+, \times} (r/\mu)$',fontsize = fs)


ax1.locator_params(axis='y', nbins = 3)
ax1.locator_params(axis='x', nbins = 5)

ax1.tick_params(axis='both', which='major', labelsize=16)







limY = 1.00
limX = 3
ax1.set_ylim(-limY,limY)
ax1.set_xlim(-limX, limX)





plt.savefig('/unsafe/tok2/PAPERS/GW_Burst/SingleWaveform.png', dpi=300)




plt.show()




