from __future__ import division
import numpy as np



#Observer location



#Particle location
Norbits = 2
r0 = 6
phi0 = np.arange(0,2*np.pi*Norbits)
u = np.linspace(0,200,len(phi0))


omega0 = np.sqrt(1/r0**3)
psi = omega0*u + phi0



m = 1
Aplus = 
expr = Aplus * np.cos(m*psi) + Bplus * np.sin(m*psi)










