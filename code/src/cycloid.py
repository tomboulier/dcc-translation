# draw an arc of cycloid

import matplotlib.pyplot as plt
import numpy as np

# parameters
R0 = 1
omega = 1
T = 2*np.pi/3
vlim = -2*R0*np.sin(omega * T/2)/T
vmin = vlim - 1
#vmin = 1
vmax = 1
Nv = 10

for v in np.linspace(vmin,vmax, num=Nv):
	# functions
	t = np.arange(-T/2, T/2, 0.01)
	s1 = -R0 * np.sin(omega * t) - (t+T/2)*v
	s2 = R0 * np.cos(omega * t)
	
	xr = R0 * np.sin(omega * T/2)
	xl = -R0 * np.sin(omega * T/2) - T * v
	y0 =  R0 * np.cos(omega * T/2)
	
	# plot figures
	plt.plot(s1, s2, '-')
	plt.plot(xr,y0,'^g', xl,y0,'sr', ms = 10)
	plt.ylim((0,np.abs(s2).max()))

plt.grid(True)
plt.axes().set_aspect('equal', 'datalim')
plt.show()
