from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
import os

#Set up plotting environment
fig = plt.figure(figsize=(10,10))
ax = plt.subplot2grid((1,1), (0,0))

#Set the path
path = os.environ['SCDir']
f0 = path+'data_lambda=0.00.txt'
f1 = path+'data_lambda=1.00.txt'


#Set the observer coordinates
theta_obs = np.pi/4
phi_obs = 0.0
nx = np.sin(theta_obs)*np.cos(phi_obs)
ny = np.sin(theta_obs)*np.sin(phi_obs)
nz = np.cos(theta_obs)


def load(f):
    data = np.loadtxt(f)
    
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    dr = x*nx + y*ny + z*nz
    dr = dr - dr[0]

    return t, dr




#Load the data
t1, dr1 = load(f0)
t2, dr2 = load(f1)


#Interpolate so they share x
finterp = interpolate.interp1d(t2,dr2,fill_value = 'extrapolate')
dr2_new = finterp(t1)


#Get the difference
R = dr2_new - dr1


ax.plot(t1,R)
ax.set_xlabel('t')
ax.set_ylabel('R')


plt.show()


