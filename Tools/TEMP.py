from __future__ import division
import matplotlib.pyplot as plt
import numpy as np








c = 3e8
T = 0.1*365*24*3600
G = 6.67e-11
MSolar = 2e30
M = 4.3e6
mu = G*M*MSolar

a = ((T**2 * mu)/(4*np.pi**2))**(1/3)
ecc = 0.9
b = a*np.sqrt(1-ecc**2)
gamma = 0 * np.pi/180




plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#Set up plotting environment
fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1))

BigE = np.linspace(0,2*np.pi,1000)
x = a*(np.cos(BigE) - ecc)
y = b*np.sin(BigE)







xp = np.cos(gamma)*x - np.sin(gamma)*y
yp = np.sin(gamma)*x + np.cos(gamma)*y


x0 = xp[0]
y0 = yp[0]





dx = xp - x0
tPH = dx/c




tNS = T*(BigE - ecc*np.sin(BigE))/(2*np.pi)





ax2.plot(BigE, tNS)
ax2.plot(BigE, tPH)














ax1.plot(x,y)
ax1.plot(xp,yp)
ax1.scatter(xp[0],yp[0])
ax1.scatter(0,0)

plt.show()










