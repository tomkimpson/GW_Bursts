from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from nfft import nfft






d = 3
col = 'C0'
fs = 20




#set up plotting environment


#Define 3d/2d plotting




if (d == 3):
    fig = plt.figure(figsize=(10,10))
    ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')
elif  (d == 2):
    fig = plt.figure(figsize=(20,10))
    ax1 = plt.subplot2grid((1,2), (0,0))
    ax2 = plt.subplot2grid((1,2), (0,1))



#PLOT MPD TRAJECTORY

MPDfile = '/unsafe/tok2/GravWaves/Trajectory.txt'


data1 = np.loadtxt(MPDfile)
x = data1[:,0]
y = data1[:,1]
z = data1[:,2]


a = np.array([x[0], y[0],z[0]])
b = np.array([x[2], y[2],z[2]])


Lvector = np.cross(a,b)
Lvector = 1000*Lvector/np.linalg.norm(Lvector)





if (d == 3):
    ax1.plot(x,y,z)
    ax1.scatter(0,0,0,c='r')
    ax1.scatter(x[0],y[0],z[0],c='C3')
    ax1.set_zlim(-1000,1000)
    ax1.set_xlabel(r'$ x [\rm r_g]$',fontsize = fs)
    ax1.set_ylabel(r'$ y [\rm r_g]$',fontsize = fs)
    ax1.set_zlabel(r'$ z [\rm r_g]$',fontsize = fs)
    
    '''
    #Plot observer vectors
    xx = [0,1000]
    yy = [0,0]
    zz = [0,1000]
    ax1.plot(xx,yy,zz)
    
    #Plot vector nomal to orbital plane
    xx = [x[0],Lvector[0]]
    yy = [y[0],Lvector[1]]
    zz = [z[0],Lvector[2]]
    ax1.plot(xx,yy,zz)
    '''
elif (d==2):
    ax1.plot(x,y,c='C0')
    ax1.scatter(x[0],y[0],c='C4')
    ax2.plot(x,z,c='C1')
    ax1.set_xlabel(r'$ x [\rm r_g]$',fontsize = fs)
    ax2.set_xlabel(r'$ x [\rm r_g]$',fontsize = fs)
    ax1.set_ylabel(r'$ y [\rm r_g]$',fontsize = fs)
    ax2.set_ylabel(r'$ z [\rm r_g]$',fontsize = fs)
    ax1.scatter(0,0,c='r')
plt.show()








plt.show()




