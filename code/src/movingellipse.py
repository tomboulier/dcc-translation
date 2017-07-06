#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
 
"""
    This script simulate a moving ellipse, illuminated by a source following
    an arc of circle. The illumination is fanbeam.
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

class Simulator(object):
	"""
        In this class, all simulations are performed
    """

	def __init__(self, params):
		"""
			Load the parameters
		"""
		self.params = params

	def run(self):
		"""
			This is the core of the computation
		"""

		# get general parameters
		R0 = self.params.R0
		sdd = self.params.sdd
		omega = self.params.omega
		T = self.params.T
		imageSize = self.params.imageSize
		Nt = self.params.Nt
		v = self.params.v
		y0 = self.params.y0

		# get ellipse parameters
		ellipseDensity = self.params.ellipseDensity
		ellipseAngle = self.params.ellipseAngle
		ellipseCenterX = self.params.ellipseCenterX
		ellipseCenterY = self.params.ellipseCenterY
		ellipseSemiAxis = [self.params.ellipseSemiAxisX,
						   self.params.ellipseSemiAxisY,
						   self.params.ellipseSemiAxisY]

		projarray = np.zeros((Nt,imageSize)) # the array giving projections

		for nt in range(Nt):
			t = -T/2 + T * nt/(Nt-1)
			angle = omega*t + 90
		
			# create geometry
			geometry = srtk.ThreeDCircularProjectionGeometry()
			geometry.AddProjection(R0, sdd, angle, 0, 0)
		
			## simulate fan-beam acquisition
			# create empty stack
			proj = srtk.ConstantImageSource()
			proj.SetSize([imageSize,1,1])
			proj.SetSpacing([1, 1, 1])
			orig = (np.array(proj.GetSize())-1)*np.array(proj.GetSpacing())*-.5
			proj.SetOrigin(orig)
			proj.SetConstant(0.0)
			source = proj.Execute()
		
			# compute intersection with ellipse
			rei = srtk.RayEllipsoidIntersectionImageFilter()
			rei.SetDensity(ellipseDensity)
			rei.SetAngle(ellipseAngle)
			rei.SetCenter([ellipseCenterX + (t+T/2)*v, ellipseCenterY, 0])
			rei.SetAxis(ellipseSemiAxis)
			rei.SetGeometry(geometry)
			reiImage = rei.Execute(source)
		
			# put this projection into projarray
			projarray[nt,:] = srtk.GetArrayFromImage(reiImage)[0,0,:]

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



