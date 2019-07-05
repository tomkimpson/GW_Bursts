import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema
import scipy
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#SEt up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))
#Load data

datafile = '/home/tok2/GW_Bursts/src/SNR_Data.txt'
data1 = np.loadtxt(datafile)
r = data1[:,0]
S = data1[:,1]

ax1.scatter(r,S)
ax1.set_yscale('log')
plt.show()

