#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
 
"""
    This script simulate a moving ellipse, illuminated by a source
    following an arc of circle. The illumination is fanbeam.
"""

import SimpleRTK as srtk
import numpy as np
import ConfigParser
import ipdb
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class Parameters(object):
	"""
        This class stores all the parameters used in the simulation
    """
 
	def __init__(self, configFile):
		"""
    		The parameters are:
    		- R0: radius of the circle, in mm
    		- sdd: source to detector distance, in mm
    		- T: speed of rotation, in degree/s
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
		self.y0 = self.R0*np.cos(self.omega*self.T/2 * 2*np.pi/360)

	def get_angle(self, t):
		return self.omega * t + 90

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
		self.params = params

	def get_geometry(self,t):
		geometry = srtk.ThreeDCircularProjectionGeometry()
		geometry.AddProjection(self.params.R0,self.params.sdd,
							   self.params.get_angle(t), 0, 0)
		return geometry

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

		for nt in range(Nt):
			t = -T/2 + T * nt/(Nt-1)
			projarray[nt,:] = self.ellipse.compute_projection(t, self.source, self.detector)

		# store results
		results = Results(self.params)
		results.projections = projarray
		return results



class Results(object):
	"""
		Encapsulation of everything that is computed by Simulator
	"""
	def __init__(self, params):
		self.params = params
		self.projections = None

	def plotSinogram(self):
		plt.imshow(self.projections[:,:], cmap = cm.Greys_r)
		plt.show()


if __name__ == '__main__':
	p = Parameters('test.ini')
	s = Simulator(p)
	res = s.run()

	res.plotSinogram()



