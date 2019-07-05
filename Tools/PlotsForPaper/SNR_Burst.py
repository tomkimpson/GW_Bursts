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
ConvertDays = 60*60*24
tplot = t / (ConvertDays) # days

if BigTheta == 0:
    hplus1 = data1[:,6]
    hcross1 = data1[:,7]
elif BigTheta == np.pi/4:
    hplus1 = data1[:,8]
    hcross1 = data1[:,9]
elif BigTheta == np.pi/2:
    hplus1 = data1[:,10]
    hcross1 = data1[:,11]






#Normalise to t=0 and only get the data within Tobs/2 of periapsis

r = data1[:,12]
ind = np.argmin(r)
rmin = r[ind]
tmin = tplot[ind]
tplot = tplot- tmin

ax1.plot(tplot,hplus1)
ax1.plot(tplot,hcross1)


Tobs = 0.500*24*60*60 / ConvertDays # 1 day observation

t_upper = tmin + Tobs/2
t_lower = tmin - Tobs/2


ax1.axvline(t_upper - tmin, linestyle = '--', c='0.5')
ax1.axvline(t_lower - tmin, linestyle = '--', c='0.5')








#Get the data witin the bounds to use for SNR calculations

tmin = t[ind]
Tobs = Tobs * ConvertDays
t_upper = tmin + Tobs/2
t_lower = tmin - Tobs/2

ConversionFactor = data1[0,15]
hplus = []
hcross = []
tt = []

for i in range(len(t)):
    if t_lower<t[i]<t_upper:
        tt.extend([t[i]])
        hplus.extend([hplus1[i]])
        hcross.extend([hcross1[i]])



# Can plot just as a sanity check
#ax1.scatter(np.array(tt)/ConvertDays - (tmin/ConvertDays),hplus)


hplus = np.array(hplus) / ConversionFactor
hcross = np.array(hplus) / ConversionFactor
t = tt




#Interpolate to get even sampling
fs = 2.0 #sampling frequency
dt = 1.0/fs
t1 = np.arange(tmin-Tobs/2, tmin+Tobs/2, dt)
hplus = np.interp(t1,t,hplus)
hcross= np.interp(t1,t,hcross)



#Get the frequencies
f = np.fft.rfftfreq(hplus.size, dt)
df = f[1] - f[0] #arethe frequencies evelyspaces?
print ('LENGTHS:', len(hplus))






#Calculate the FT
hplusT = dt*np.fft.rfft(hplus) #/ factorW
hcrossT = dt*np.fft.rfft(hcross) #/ factorW






print ('Got the FT')






#Get rid of zeroth frequencies

hplusT = hplusT[1:] # get rid of zeroth frequency
hcrossT = hcrossT[1:]
f = f[1:]
print ('Completed FT')



#LISA NOISE CURVE
Larm = 2.5e9
Clight = 3e8
fstar = Clight/(2*np.pi*Larm)
NC = 2

#Constants for glactic binary confusion noise
alpha = 0.133
beta = 243.
kappa = 482.
gamma = 917.
f_knee = 2.58e-3

A = 1.8e-44/NC
Sc = 1. + np.tanh(gamma*(f_knee-f))
Sc *=np.exp(-f**alpha + beta*f*np.sin(kappa*f))

Sc *= A*f**(-7./3.)


#Response function load from folder https://github.com/eXtremeGravityInstitute/LISA_Sensitivity 
RFILE = np.loadtxt('/unsafe/tok2/ResponseFunction.txt')
Rx = RFILE[:,0] * fstar
Ry = RFILE[:,1] * NC

newR = np.interp(f,Rx,Ry)
R = newR




#R = 0.30/(1+0.6*(f/fstar)**2)


#Power Spectral density
P_oms = (1.5e-11)**2 * (1. + (2.0e-3/f)**4)
P_acc = (3.0e-15)**2 * (1.+(0.4e-3/f)**2)*(1. + (f/(8.0e-3))**4)
Pn = (P_oms + 2.*(1. + np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4)/Larm**2



#Total noise
S = Pn/R + Sc

#S = 10/(3*Larm**2) * (P_oms + 2*(1+np.cos(f/fstar)**2)*P_acc/(2*np.pi*f)**4) * (1+0.60*(f/fstar)**2)



hsig = abs(hplusT)**2 + abs(hcrossT)**2


SNR2 = 4 * scipy.integrate.simps((hsig)/S , f)
SNR = np.sqrt(SNR2)
print ('The calculated SNR = ', SNR)





#ax2.set_ylim(1e-22,1.0e-13)


#Make it look pretty

fs = 20
#AX1
ax1.set_xlabel(r' t [days]',fontsize = fs)
ax1.set_ylabel(r'$ h_{+, \times} (r/\mu)$',fontsize = fs)
ax1.locator_params(axis='y', nbins = 3)
ax1.locator_params(axis='x', nbins = 5)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_xlim(-4,4)
#ax1.set_ylim(-0.7,0.7)

#AX2
ax2.loglog(f,np.sqrt(hsig), C = 'C3')
ax2.loglog(f,np.sqrt(S), C='C2')
ax2.set_xlim(1.0e-5,1.0e0)
ax2.set_xlabel(r' $f$ [Hz]',fontsize = fs)
ax2.set_ylabel(r'$ \tilde{h}(f)$ [Hz] $^{-1}$',fontsize = fs)
ax2.tick_params(axis='both', which='major', labelsize=16)




#Full figure save
plt.savefig('/unsafe/tok2/PAPERS/GW_Burst/SNR_TOT.png', dpi=300)


bbox = ax1.get_tightbbox(fig.canvas.get_renderer())
extent = bbox.transformed(fig.dpi_scale_trans.inverted())
plt.savefig("/unsafe/tok2/PAPERS/GW_Burst/SNRa_47TUC_INC.png".format(ax1),bbox_inches=extent.expanded(1.1,1))



bbox = ax2.get_tightbbox(fig.canvas.get_renderer())
extent = bbox.transformed(fig.dpi_scale_trans.inverted())
plt.savefig("/unsafe/tok2/PAPERS/GW_Burst/SNRb_47TUC_INC.png".format(ax2),bbox_inches=extent.expanded(1.1,1))

plt.show()




