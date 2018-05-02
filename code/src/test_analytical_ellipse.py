#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
 
"""
    This script computes analytically the projections of
    a moving ellipse, illuminated by a source which follows
    an arc of circle. The illumination is fanbeam.
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
    		- imageSize: length of the detector, in mm
    		- y0: distance from isocenter to the chord of the
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
		self.v2 = cfg.getfloat('Parameters', 'v2')
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

class AnalyticalEllipse(object):
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
		v2 = self.params.v2

		return [self.params.ellipseCenterX + (t+T/2)*v,
		        self.params.ellipseCenterY + (t+T/2)*v2,
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
		N = self.params.imageSize
		results = np.zeros(N)

		# general parameters
		alpha = self.params.get_alpha_range()
		omega = self.params.omega/360 * 2*np.pi

		# ellipse parameters
		x,y,_ = self.get_center(t)
		a,b,_ = self.get_axis()
		if (a != b):
			raise ValueError("Ellipse is not a circle (the analytical formula only works with circle)", a, b)
		
		s1,s2 = source.get_position(t)

		# for i in np.arange(N):
		# 	# i is the number of the pixel in the image printed on the detector
		# 	# there is no "resolution" parameter, meaning that there is 1 pixel
		# 	# per millimeter
		# 	# TODO : add this parameter?
			
		# 	# phi is the angle between the beam and the y-axis
		# 	phi = omega*t + alpha[i]

		# 	# computation of the distance between the center of the circle
		# 	# and the beam
		# 	dist = np.abs( (s1-x)*np.cos(phi) + (s2-y)*np.sin(phi) )

		# 	# stores in the array
		# 	if dist > a:
		# 		results[i] = 0
		# 	else:
		# 		results[i] = 2 * np.sqrt(a**2 - dist**2)

		# phi is the angle between the beam and the y-axis
		phi = omega*t + alpha
		
		# distance between the center of the circle and the beam
		dist = np.abs( (s1-x)*np.cos(phi) + (s2-y)*np.sin(phi) )

		# results = (dist<a) * 2 * np.sqrt(a**2 - dist**2)
		# [x if x < 5 else 0 for x in np.arange(10)]
		# ipdb.set_trace()
		results[dist < a] = (2 * np.sqrt(a**2 - dist**2))[dist<a]

		

		return self.get_density() * results

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
		self.ellipse = AnalyticalEllipse(params)

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

if __name__ == '__main__':

	p = Parameters('circle_test.ini')
	res = Simulator(p).run()
	res.plotSinogram()

















