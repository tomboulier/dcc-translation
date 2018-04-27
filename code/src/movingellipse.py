#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
 
"""
    This script simulate a moving ellipse, illuminated by a source
    following an arc of circle. The illumination is fanbeam.
"""
import ipdb

import math
import SimpleRTK as srtk
import numpy as np
import ConfigParser
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
matplotlib.rcParams['text.latex.unicode'] = True
import matplotlib.cm as cm
from scipy import interpolate
from scipy.spatial import distance
from scipy import integrate
from sklearn.metrics import mean_squared_error
from math import sqrt

class Parameters(object):
	"""
        This class stores all the parameters used in the simulation
    """
 
	def __init__(self, configFile):
		"""
    		The parameters are:
    		- R0: radius of the circle, in mm
    		- sdd: source to detector distance, in mm
    		- omega: speed of rotation, in degree/s
    		- T: duration of acquisition, in s
    		- y0: total time of experiment, in s
    		- imageSize: distance from isocenter to the chord of the
    		             arc of circle, in mm
    		- Nt: number of discretization points for the time
    		- v: speed of the object along the x-axis, in mm/s
    	"""
		cfg = ConfigParser.ConfigParser()
		cfg.read(configFile)
		# general parameters
		self.R0 = cfg.getfloat('Parameters', 'R0')
		self.sdd = cfg.getfloat('Parameters', 'sdd')
		self.omega = cfg.getfloat('Parameters', 'omega')
		self.T = cfg.getfloat('Parameters', 'T')
		self.imageSize = cfg.getint('Parameters', 'imageSize')
		self.Nt = cfg.getint('Parameters', 'Nt')
		self.v = cfg.getfloat('Parameters', 'v')
		# parameters describing the ellipse
		self.ellipseDensity = cfg.getfloat('Ellipse', 'density')
		self.ellipseAngle = cfg.getfloat('Ellipse', 'angle')
		self.ellipseCenterX = cfg.getfloat('Ellipse', 'x')
		self.ellipseCenterY = cfg.getfloat('Ellipse', 'y')
		self.ellipseSemiAxisX = cfg.getfloat('Ellipse', 'a')
		self.ellipseSemiAxisY = cfg.getfloat('Ellipse', 'b')

		# These parameters are derived from previous ones
		self.y0 = self.R0*np.cos(self.omega*self.T/2 / 360 *2*np.pi)

	def get_angle(self, t, unit = 'degrees'):
		return self.omega * t + 90

	def get_max_angle(self):
		return self.get_angle(self.T/2)

	def get_min_angle(self):
		return self.get_angle(-self.T/2)

	def get_time_range(self):
		nt = np.arange(self.Nt)
		return -self.T/2 + self.T * nt/(self.Nt-1)

	def get_alpha_range(self):
		x = np.arange(-self.imageSize/2, self.imageSize/2)
		# return np.arctan( x / self.sdd) * 360 / (2*np.pi)
		return np.arctan( x / self.sdd)

	def get_virtual_positions_vector(self):
		"""
			These will be the abciss of points where DCC function
			will be computed.
		"""
		xmax = self.R0 * np.sin(self.omega*self.T/2 /360*2*np.pi) * .75
		return np.linspace(-xmax, xmax)

class MovingEllipse(object):
	"""docstring for MovingEllipse"""
	def __init__(self, params):
		self.params = params

	def get_density(self):
		return self.params.ellipseDensity

	def get_angle(self):
		return self.params.ellipseAngle

	def get_center(self, t):
		"""
			Since it is moving, the position depends on t
		"""
		T = self.params.T
		v = self.params.v

		return [self.params.ellipseCenterX + (t+T/2)*v,
		        self.params.ellipseCenterY,
		        0]

	def get_axis(self):
		return [self.params.ellipseSemiAxisX,
				self.params.ellipseSemiAxisY,
				self.params.ellipseSemiAxisY]

	def compute_projection(self,t,source,detector):
		"""
			Simulate fan-beam acquisition of the object with given
			source and detector, at time t
		"""
		# create geometry of the source at time t
		geometry = source.get_geometry(t)

		# compute intersection of fan-beam with ellipse
		empty_image_detector = detector.get_empty_image()
		rei = srtk.RayEllipsoidIntersectionImageFilter()
		rei.SetDensity(self.get_density())
		rei.SetAngle(self.get_angle())
		rei.SetCenter(self.get_center(t)) # 
		rei.SetAxis(self.get_axis())
		
		rei.SetGeometry(geometry)
		reiImage = rei.Execute(empty_image_detector)
	
		return srtk.GetArrayFromImage(reiImage)[0,0,:]

class Detector(object):
	"""docstring for Detector"""
	def __init__(self, params):
		self.params = params

	def get_empty_image(self):
		proj = srtk.ConstantImageSource()
		proj.SetSize([self.params.imageSize,1,1])
		proj.SetSpacing([1, 1, 1])
		origin = (np.array(proj.GetSize()) - 1) \
		         * np.array(proj.GetSpacing()) \
		         * -.5
		proj.SetOrigin(origin)
		proj.SetConstant(0.0)
		return proj.Execute()

class Source(object):
	"""docstring for Source"""
	def __init__(self, params):
		self.R0 = params.R0
		self.sdd = params.sdd
		self.omega = params.omega
		self.get_angle = lambda t: params.get_angle(t)

	def get_geometry(self,t):
		geometry = srtk.ThreeDCircularProjectionGeometry()
		geometry.AddProjection(self.R0,self.sdd,
							   self.get_angle(t), 0, 0)
		return geometry

	def get_position(self, t):
		return np.array((-self.R0 * np.sin(self.omega * t / 360 * 2*np.pi),\
			              self.R0 * np.cos(self.omega * t / 360 * 2*np.pi)))

class Simulator(object):
	"""
        In this class, all simulations are performed
    """

	def __init__(self, params):
		"""
			Load the parameters
		"""
		self.params = params
		self.source = Source(params)
		self.detector = Detector(params)
		self.ellipse = MovingEllipse(params)

	def run(self):
		"""
			Computes all the projections, at any angle and any time
		"""
		# get general parameters
		T = self.params.T
		imageSize = self.params.imageSize
		Nt = self.params.Nt

		# the array giving projections
		projarray = np.zeros((Nt,imageSize))

		for nt,t in enumerate(self.params.get_time_range()):
			projarray[nt,:] = self.ellipse.compute_projection(t,
				                                              self.source,
				                                              self.detector)

		# store results
		results = Results(self.params, self.source, self.detector)
		results.projections = projarray
		return results

class Results(object):
	"""
		Encapsulation of everything that is computed by Simulator
	"""
	def __init__(self, params, source, detector):
		self.params = params
		self.source = source
		self.detector = detector
		self.projections = None
		self.projections_interpolator = None
		self.DCC_function = None
		self.DCC_function_theo = None

	def plotSinogram(self, xunits = 'mm'):
		# define the limits of the axis
		imageSize = self.params.imageSize
		T = self.params.T
		max_angle = self.params.get_max_angle()
		min_angle = self.params.get_min_angle()
		phimax = np.arctan( .5*imageSize / self.params.sdd) * 360 / (2*np.pi)

		# plot the image
		plt.figure()
		
		if xunits == 'mm':
			# the units here represent a distance (on the detector)
			plt.xlabel('Distance from detector center (in mm)', labelpad=20)
			extent = [-imageSize/2, imageSize/2, max_angle, min_angle]
			aspect = imageSize / (max_angle - min_angle)

		elif xunits == 'degrees':
			# the units here represent an angle ('phi' in T(x,phi))
			plt.xlabel('Beam direction (in degrees)', labelpad=20)
			extent = [-phimax, phimax, max_angle, min_angle]
			aspect = 2 * phimax / (max_angle - min_angle)

		plt.imshow(self.projections, cmap = cm.Greys_r, extent = extent, 
			       aspect = aspect/2)
		plt.ylabel('Gantry angle (in degrees)', labelpad=20)
		matplotlib.rcParams.update({'font.size': 22})
		plt.show()

	def interpolate_projection(self):
		"""
			Interpolation of the operator T(alpha,t).
			
			Be careful: the angle alpha is the angle between the beam and the
			line joining the source to the center of the detector. Not to be
			confused with phi, which is the angle between the beam and the 
			y-axis
		"""
		t = self.params.get_time_range()
		alpha = self.params.get_alpha_range()

		self.projections_interpolator = interpolate.interp2d(alpha, t,
			                                                 self.projections,
			                                                 kind='linear')

class DataConsistencyConditions(object):
	"""
		The computations of DCCs happen here
	"""
	def __init__(self, results):
		self.results = results
		self.params = results.params
		self.source = results.source

		self.DCC_function = None

	def get_virtual_source_position(self, t, v):
		"""
			Computes :math:`s_v(t)`
		"""
		Mvt = np.array( ( (t+self.params.T/2) * v, 0) )
		return self.source.get_position(t) - Mvt

	def alpha(self, t, x, v):
		"""
			Since 'lambda' is a reserved keword in Python, we use 'alpha'
			to express the function :math:`\lambda` in Theorem 1, i.e.
			
			.. math::
			\lambda_n(t,x) = \arctan \left( \frac{x + R_0 \sin(\omega t) + \
			\left( t + \frac{T}{2} \right)v}{R_0 \cos(\omega t) - y_0} \
			\right)
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		y0 = self.params.y0
	
		num = x + R0 * np.sin(omega*t) + (t+T/2) * v
		denom = R0 * np.cos(omega*t) - y0
		
		if denom==0.0:
			return np.sign(num) * np.pi / 2
		else:
			return np.arctan(num/denom)

	def jacobian(self, t, x, v):
		"""
			The Jacobian used in the change of variable in Theorem 1,
			which is given by

			..math::
			J(x,t,v) = \frac{ R_0^2 \omega - v y_0 + R_0 \cos(\omega \
			t)(v-\omega y_0) + R_0 \omega \sin(\omega t)(x + \left( t \
			+ \frac{T}{2} \right)v ) }{ \left( R_0 \cos(\omega t) - \
			y_0 \right)^2 } dt
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		y0 = self.params.y0

		num = R0**2*omega - v*y0 + R0 * np.cos(omega*t) * (v - omega*y0) \
		                         + R0 * omega * np.sin(omega*t) * (x + (t+T/2)*v)
		denom = (R0 * np.cos(omega*t) - y0)**2

		return num / denom

	def W(self, t, x, n, v):
		"""
			The weighting function in Theorem 1, i.e.

			..math::
			W_n(x,t,v) = \frac{ \left( x+R_0 \sin(\omega t) + \left( \
			t + \frac{T}{2} \right)v \right)^n }{D_{x,t} \left( R_0 \
			\cos(\omega t) - y_0 \right)^{n-1}} J(x,t,v) dt,

			where :math:`J(x,t,v)` is given by the function 'jacobian'
		"""
		R0 = self.params.R0
		omega = self.params.omega
		T = self.params.T
		y0 = self.params.y0

		s_v_t = self.get_virtual_source_position(t, v)
		D_x_t = distance.euclidean(s_v_t, (x,y0))

		num = (x + R0 * np.sin(omega*t/360 *2*np.pi) + (t+T/2)*v)**n
		denom = D_x_t * (R0 * np.cos(omega*t/360 *2*np.pi) - y0)**(n-1)

		return num / denom * self.jacobian(t,x,v)

	def integrand(self, t, x, n, v):
		"""
			The function to be integrated in DCC, i.e.

			..math::
			T(t,\lambda_t) W_n(t,x)

			We put it in a function since we have to be careful with Nans
		"""

		T = self.params.T
		omega = self.params.omega / 360 * 2*np.pi

		alpha_v = lambda t,x: self.alpha(t, x, v)
		fb_proj = lambda t,x: self.results.projections_interpolator(alpha_v(t,x)-omega*t, t)
		weight = lambda t,x,n: self.W(t, x, n, v)

		y = fb_proj(t,x) * weight(t,x,n)

		if fb_proj(t,x) < 1e-10:
			return 0.
		
		if math.isnan(y):
			return 0.

		return y

	def B(self, x, n, v, tol = 1):
		"""
			Compute the function B_n(x) in Theorem 1, which is supposed
			to be a polynom of order at most n, where n is the order of DCC.
		"""
		T = self.params.T
		Nt = self.params.Nt

		integrand = lambda t,x,n: self.integrand(t,x,n,v)

		t = np.linspace(-T/2, T/2, Nt)
		y = np.array([integrand(time,x,n) for time in t])
		
		return integrate.simps(y, dx = t[1]-t[0])

	def compute_DCC_function(self, v, n):
		"""
			Transform B as a lambda function
		"""
		if self.results.projections_interpolator is None:
			self.results.interpolate_projection()

		self.DCC_function = lambda x,n: self.B(x, n, v)

	def compute_vectorized_DCC_function(self, v, x, n):
		"""
			Compute a vector giving all values of B(x) for each
			point in x
		"""
		self.compute_DCC_function(v,n)
		Bn = np.vectorize(lambda x: self.DCC_function(x, n))

		return Bn(x)

class Optimizator(object):
	"""
		In order to avoid inverse crime, we have to be
		very careful here with the parameters the DCCs
		can access.
	"""

	def __init__(self, results, dcc):
		self.results = results
		self.dcc = dcc

	def residual_polyfit(self, x, n, v):
		"""
			Computes the function |Bn(x)-P(x)|, where P(x) is the
			polynom obtained by fitting.
		"""

		# self.dcc.compute_DCC_function(v)
		# Bn = np.vectorize(lambda x: self.dcc.DCC_function(x, n))
		Bn = self.dcc.compute_vectorized_DCC_function(v,x,n)
		_, res, _, _, _ = np.polyfit(Bn(x), y, n, full = True)

		return res[0]	

	
if __name__ == '__main__':

	p = Parameters('example.ini')
	s = Simulator(p)
	res = s.run()
	# res.plotSinogram()

	# Plot DCC
	v = p.v
	n = 2
	x = p.get_virtual_positions_vector()

	# compute x -> Bn(x) function
	DCC = DataConsistencyConditions(res)
	y = DCC.compute_vectorized_DCC_function(v,x,n)

	# interpolation with polynom
	poly = np.polyfit(x, y, n)
	yfit = np.poly1d(poly)(x)
	# rmse = sqrt(mean_squared_error(y, yfit))
	diff = y-yfit
	rmse = np.sqrt((diff*diff).sum())/np.sqrt((y*y).sum())
	textrmse = r"$%.4f$" % (rmse)
	textstr = r"$RMSE = $" + textrmse

	# plot results
	fig, ax = plt.subplots(1)
	plt.plot(x, y, 'ob')
	plt.plot(x, yfit, '-r')
	axes = plt.gca()
	axes.set_ylim([y.mean()-20,y.mean()+20])
	matplotlib.rcParams.update({'font.size': 25})
	ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30,
        verticalalignment='top')
	plt.savefig('B' + str(n) + '.eps')
	plt.show()

	# # optimization
	# print "Error of interpolation is: " + str(res.residual_polyfit(x,n,v))
	# from scipy.optimize import minimize
	# # residual_callable = lambda v: res.residual_polyfit(x,0,v) + res.residual_polyfit(x,1,v) + res.residual_polyfit(x,2,v) + res.residual_polyfit(x,3,v)
	# residual_callable = lambda v: res.residual_polyfit(x,n,v)
	# min_v = minimize(residual_callable, 0, method = 'Powell')
	# print "*** Results of optimization ***"
	# print "True value is: " + str(p.v)
	# print "Search result is: " + str(min_v.x)
	# print "Difference is: " + str(min_v.x-v)
