#!/usr/bin/env python
# -*- coding: utf-8 -*-
 
 
"""
    This script simulates a moving ellipse, illuminated by a source
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

	def get_virtual_source_position(self, t, v1, v2):
		"""
			Computes :math:`s_v(t)` given by

			..math::
			s_{\mathbf{v}}(t) = \mathcal{R}_{-\beta} \left( s(t) - M_v(t) \right)
		"""
		# translation vector
		Mvt = (t+self.params.T/2) * np.array( (v1, v2) )

		# rotation matrix
		beta = self.beta(v1,v2)
		c, s = np.cos(-beta), np.sin(-beta)
		R_beta = np.array(((c,-s), (s, c)))
		
		return R_beta.dot(self.source.get_position(t) - Mvt)

	def get_virtual_points_vector(self, v1, v2):
		"""
			These will be the x-coordinates of points where DCC function
			will be computed.
		"""
		# xmax = self.R0 * np.sin(self.omega*self.T/2 /360*2*np.pi) * .75
		# return np.linspace(-xmax, xmax)

		# take notations of Definition 3 for x', i.e. x' is between x'_l and x'_r
		# whith x'_l = s_v(T/2) and x'_r = s_v(-T/2)
		T = self.params.T
		# v1 = self.params.v
		# v2 = self.params.v2

		xl = self.get_virtual_source_position(T/2, v1, v2)[0]
		xr = self.get_virtual_source_position(-T/2, v1, v2)[0]

		return np.linspace(xl,xr)

	def beta(self, v1, v2):
		"""
			The rotation angle :math:`\beta` defined in formula (7)
			of the abstract

			.. math::
			\beta = \arctan \left( \frac{T v_2} \
			{2 R_0 \sin(\omega T/2) + T v_1} \right)
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T

		num = T * v2
		denom = 2 * R0 * np.sin(omega*T/2) + T*v1

		if denom == 0.:
			raise ZeroDivisionError("denominator is equal to zero in function beta")
		else:
			return np.arctan(num/denom)

	def lambda_v(self, t, v1, v2, x):
		"""
			Function :math:`\lambda^{(v)}` in Theorem 1, i.e.
			
			..math::
			\lambda^{(\bv)}(t,x') = \arctan \left( F^{(\bv)}(t,x') \right),
			
			where :math:`F^{(\bv)}(t,x')` is defined as the fraction
			:math:`A/B` with
			
			..math::
				A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) + \
				\sin \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right)
			
			and
			..math::
				B = \cos \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right) - \
				\sin \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) - y'_0 
			
		"""
		
		return np.arctan(self.F(t,v1,v2,x))

	def F(self, t, v1, v2, x):
		"""
			Function :math:`F^{(\bv)}(t,x')` which is defined in
			Theorem 1 as the fraction :math:`A/B` with
			
			..math::
				A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) + \
				\sin \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right)
			
			and
			..math::
				B = \cos \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right) - \
				\sin \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) - y'_0 
			
		"""

		A = self.F_numerator(t, v1, v2, x)
		B = self.F_denominator(t, v1, v2, x)
		
		if B ==0.:
			return np.sign(A) * np.inf
		else:
			return A/B

	def F_numerator(self, t, v1, v2, x):
		"""
			The numerator of the function :math:`F^{(\bv)}(t,x')`
			defined in Theorem, where it is denoted :math:`A`, with

			..math::
				A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) + \
				\sin \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right)
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		beta = self.beta(v1,v2)

		return x + np.cos(beta) * (R0 * np.sin(omega*t) + (t+T/2)*v1) + \
			np.sin(beta) * (R0 * np.cos(omega*t) - (t+T/2)*v2)

	def F_denominator(self, t, v1, v2, x):
		"""
			The denominator of the function :math:`F^{(\bv)}(t,x')`
			defined in Theorem, where it is denoted :math:`B`, with

			..math::
				B = \cos \beta \left( R_0 \cos(\omega t) - \
				\left( t + T/2 \right)v_2 \right) - \
				\sin \beta \left( R_0 \sin(\omega t) + \
				\left( t + T/2 \right)v_1 \right) - y'_0 
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		beta = self.beta(v1,v2)
		y0_prime = R0 * np.cos(omega*T/2 + beta)

		return np.cos(beta) * (R0 * np.cos(omega*t) - (t+T/2)*v2) - \
			np.sin(beta) * (R0 * np.sin(omega*t) + (t+T/2)*v1) - y0_prime

	def diff_t_F_numerator(self, t, v1, v2, x):
		"""
			Derivative of the function :math:`A` defined in equation (13),
			with respect to time.
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		beta = self.beta(v1,v2)

		return np.cos(beta) * (R0 * omega * np.cos(omega*t) + v1) + \
			np.sin(beta) * (- R0 * omega * np.sin(omega*t) - v2)

	def diff_t_F_denominator(self, t, v1, v2, x):
		"""
			Derivative of the function :math:`B` defined in equation (14),
			with respect to time.
		"""
		R0 = self.params.R0
		omega = self.params.omega/360 *2*np.pi
		T = self.params.T
		beta = self.beta(v1,v2)

		return np.cos(beta) * (- R0 * omega * np.sin(omega*t) - v2) - \
			np.sin(beta) * (R0 * omega * np.cos(omega*t) + v1)

	def jacobian(self, t, v1, v2, x):
		"""
			The Jacobian used in the change of variable in Theorem 1,
			which is given by the derivative of :math:`\lambda^{(\mathbf{v})}(t,x')`
			with respect to :math:`t`. Hence, one has

			..math::
			\partial_t \lambda^{(\mathbf{v})} = \
			\frac{\partial_t A B - A \partial_t}{B^2} \frac{1}{1 + F^2}
		"""
		A = self.F_numerator(t, v1, v2, x)
		B = self.F_denominator(t, v1, v2, x)
		diff_t_A = self.diff_t_F_numerator(t, v1, v2, x)
		diff_t_B = self.diff_t_F_denominator(t, v1, v2, x)
		F = self.F(t, v1, v2, x)

		# return (diff_t_A * B - A * diff_t_B) / (B*B * (1+F*F))
		return (diff_t_A * B - A * diff_t_B) / (B*B)

	def W(self, n, t, v1, v2, x):
		"""
			The weighting function in Theorem 1, i.e.

			..math::
			W_n^{(v)}(t,x') = \tan^n \left( \lambda^{(v)}(t,x') \right) \\
			\cos \left( \lambda^{(v)}(t,x') \right) J(x,t,v),

			where :math:`J(x,t,v)` is given by the function 'jacobian'
		"""
		lambda_v = self.lambda_v(t, v1, v2, x)
		jacobian = self.jacobian(t, v1, v2, x)

		return (np.tan(lambda_v))**n * np.cos(lambda_v) * jacobian

	def integrand(self, n, t, v1, v2, x):
		"""
			The function to be integrated in DCC, i.e.

			..math::
			\mathcal{F}(t,\lambda^{(v)}(t,x)) W_n^{(v)}(t,x)

			We put it in a function since we have to be careful with Nans
		"""

		T = self.params.T
		omega = self.params.omega / 360 * 2*np.pi
		beta = self.beta(v1,v2)

		lambda_v = lambda t,x: self.lambda_v(t,v1,v2,x)
		# fb_proj = lambda t,x: self.results.projections_interpolator(lambda_v(t,x)-omega*t, t)
		fb_proj = lambda t,x: self.results.projections_interpolator(lambda_v(t,x) + beta - omega*t, t)
		weight = lambda t,x,n: self.W(n,t,v1,v2,x)

		y = fb_proj(t,x) * weight(t,x,n)

		if fb_proj(t,x) < 1e-10:
			return 0.
		
		if math.isnan(y):
			return 0.

		return y

	def B(self, n, v1, v2, x):
		"""
			Compute the function B_n(x) in Theorem 1, which is supposed
			to be a polynom of order at most n, where n is the order of DCC.
		"""
		T = self.params.T
		Nt = self.params.Nt

		integrand = lambda n,t,x: self.integrand(n,t,v1,v2,x)

		t = np.linspace(-T/2, T/2, Nt)
		y = np.array([integrand(n,time,x) for time in t])
		
		return integrate.simps(y, dx = t[1]-t[0])

	def compute_DCC_function(self, n, v1, v2):
		"""
			Transform B as a lambda function
		"""
		if self.results.projections_interpolator is None:
			self.results.interpolate_projection()

		self.DCC_function = lambda n,x: self.B(n, v1, v2, x)

	def compute_vectorized_DCC_function(self, n, v1, v2, x):
		"""
			Compute a vector giving all values of B(x) for each
			point in x
		"""
		self.compute_DCC_function(n,v1,v2)
		Bn = np.vectorize(lambda x: self.DCC_function(n,x))

		return Bn(x)

class PolynomProjector(object):
	"""
		In this class, we fit the DCC function onto a polynomial
		function
	"""

	def __init__(self, dcc):
		self.dcc = dcc

	def fit_dcc_polynom(self, n, v1, v2, x):
		"""
			Interpolates the DCC function with a polynom of
			degree n
		"""

		# computes the array with Bn(x) values
		y = self.dcc.compute_vectorized_DCC_function(n,v1,v2,x)

		# interpolation with polynom
		poly = np.polyfit(x, y, n)


		# the results is an array of values
		return np.poly1d(poly)(x)

	def compute_RMSE(self, n, v1, v2, x):
		"""
			Compute the root mean square error of the fitting
			with fit_dcc_polynom
		"""

		y = self.dcc.compute_vectorized_DCC_function(n,v1,v2,x)
		yfit = self.fit_dcc_polynom(n,v1,v2,x)
		
		diff = y - yfit
		
		return np.sqrt((diff*diff).sum())/np.sqrt((y*y).sum())

	def plot_fitting(self, n, v1, v2, x):
		"""
			Plot the result of the fitting procedure, showing
			the function Bn(x) with the polynom and the RMSE
		"""
		# text for RMSE
		rmse = self.compute_RMSE(n,v1,v2,x)
		textrmse = r"$%.4f$" % (rmse)
		textstr = r"$RMSE = $" + textrmse

		y = self.dcc.compute_vectorized_DCC_function(n,v1,v2,x)
		yfit = self.fit_dcc_polynom(n,v1,v2,x)

		fig, ax = plt.subplots(1)
		plt.plot(x, y, 'ob')
		plt.plot(x, yfit, '-r')
		axes = plt.gca()
		# axes.set_ylim([y.mean()-20,y.mean()+20])
		matplotlib.rcParams.update({'font.size': 25})
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30,
    	    verticalalignment='top')
		# plt.savefig('B' + str(n) + '.eps')
		plt.show()

	def residual_polyfit(self, n, v1, v2, x):
		"""
			Computes the function |Bn(x)-P(x)|, where P(x) is the
			polynom obtained by fitting.
		"""
		Bn = self.dcc.compute_vectorized_DCC_function(n,v1,v2,x)
		_, res, _, _, _ = np.polyfit(x, Bn, n, full = True)

		return res[0]

class Optimizer(object):
	"""
		Optimization process in order to recover velocity 
		occurs in this class
	"""
	def __init__(self, polynom_projector):
		"""
			In order to avoid inverse crime, we have to be
			very careful here with the parameters the DCCs
			can access.
		"""
		self.polyproj = polynom_projector
		
	def minimize_rmse(self, n, x):
		"""
			Minimization of RMSE of the fitting
		"""
		from scipy.optimize import minimize

		residual_callable = lambda v: self.polyproj.residual_polyfit(n,v,x)
		# residual_callable = lambda v: self.polyproj.compute_RMSE(n,v,x)

		min_v = minimize(residual_callable, 0, method = 'Powell')

		return min_v.x
	
	
if __name__ == '__main__':

	p = Parameters('example.ini')
	s = Simulator(p)
	res = s.run()
	res.plotSinogram()

	# Plot DCC
	v = p.v
	n = 3

	DCC = DataConsistencyConditions(res)
	# DCC.compute_vectorized_DCC_function(n,p.v,p.v2,x)

	polyproj = PolynomProjector(DCC)
	x = DCC.get_virtual_points_vector(p.v,p.v2)
	# polyproj.fit_dcc_polynom(n,p.v,p.v2,x)
	# print "Error of interpolation is: " + str(polyproj.residual_polyfit(n,v,x))
	polyproj.plot_fitting(n, p.v, p.v2, x)

	# # Optimization
	# optim = Optimizer(polyproj)
	# min_v = optim.minimize_rmse(n,x)

	# print "*** Results of optimization ***"
	# print "True value is: " + str(p.v)
	# print "Search result is: " + str(min_v)
	# print "Difference is: " + str(min_v-v)
