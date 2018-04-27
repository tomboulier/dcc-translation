from movingellipse import *
from numpy.linalg import norm
import pytest

def wiki_example_RTK(configFile):
	"""
		Simulates a cone-beam illumination of an ellipse.

		Taken from RTK's wiki, see the following page:
		http://wiki.openrtk.org/index.php/SimpleRTK#Testing_SimpleRTK
	"""
	# Read the parameters from file
	cfg = ConfigParser.ConfigParser()
	cfg.read(configFile)
	numberOfProjections = cfg.getint('Parameters', 'Nt')
	omega = cfg.getfloat('Parameters', 'omega')
	T = cfg.getfloat('Parameters', 'T')
	firstAngle = 90 - omega * T/2
	angularArc = omega*T
	sid = cfg.getfloat('Parameters', 'R0')  # source to isocenter distance in mm
	sdd = cfg.getfloat('Parameters', 'sdd')  # source to detector distance in mm
	isox = 0  # X coordinate on the projection image of isocenter
	isoy = 0  # Y coordinate on the projection image of isocenter
	imageSize = cfg.getint('Parameters', 'imageSize')

	# get ellipse parameters
	ellipseDensity = cfg.getfloat('Ellipse', 'density')
	ellipseAngle = cfg.getfloat('Ellipse', 'angle')
	ellipseCenterX = cfg.getfloat('Ellipse', 'x')
	ellipseCenterY = cfg.getfloat('Ellipse', 'y')
	ellipseSemiAxisX = cfg.getfloat('Ellipse', 'a')
	ellipseSemiAxisY = cfg.getfloat('Ellipse', 'b')

	# Defines the RTK geometry object
	geometry = srtk.ThreeDCircularProjectionGeometry()
	for x in range(0, numberOfProjections):
	    angle = firstAngle + x * angularArc / (numberOfProjections-1)
	    geometry.AddProjection(sid, sdd, angle, isox, isoy)
	
	# initialization of 'blank' image
	constantImageSource = srtk.ConstantImageSource()
	sizeOutput = [ imageSize, imageSize,  numberOfProjections ]
	spacing = [ 1.0, 1.0, 1.0 ]
	origin = (np.array(sizeOutput)-1) * np.array(spacing) * -.5
	constantImageSource.SetOrigin( origin )
	constantImageSource.SetSpacing( spacing )
	constantImageSource.SetSize( sizeOutput )
	constantImageSource.SetConstant(0.0)
	source = constantImageSource.Execute()
	 
	rei = srtk.RayEllipsoidIntersectionImageFilter()
	semiprincipalaxis = [ ellipseSemiAxisX,
						  ellipseSemiAxisY,
						  ellipseSemiAxisY]
	center = [ ellipseCenterX, ellipseCenterY, 0]
	# Set GrayScale value, axes, center...
	rei.SetDensity(ellipseDensity)
	rei.SetAngle(ellipseAngle)
	rei.SetCenter(center)
	rei.SetAxis(semiprincipalaxis)
	rei.SetGeometry( geometry )
	reiImage = rei.Execute(source)

	return srtk.GetArrayFromImage(reiImage[:,imageSize/2,:])

configFile = 'test.ini'
tol = .01

# simulation by our algorithm
p = Parameters(configFile)
s = Simulator(p)
res = s.run()

# simulation with RTK's example in wiki
projWiki = wiki_example_RTK(configFile)

# difference between the two projections
diff = res.projections - projWiki

## General tests
def test_get_angle():
	"""
		The angles start at 90 since we compute the angle with the
		y-axis
	"""
	assert p.get_angle(0) == 90

def test_y0():
	"""
		:math:`y0` is supposed to be stricly positive
	"""
	assert p.y0 > 0

## Testing projection operator

def test_projections_have_proper_dimensions():
    assert res.projections.shape[0] == p.Nt
    assert res.projections.shape[1] == p.imageSize

def test_projections_fits_example_in_wiki_RTK():
	"""
		Compare the two projections, in L2 norm
	"""
	assert norm(diff) <= tol * norm(projWiki)
	assert norm(diff,1) <= tol * norm(projWiki,1)
	assert norm(diff,2) <= tol * norm(projWiki,2)
	assert norm(diff,np.inf) <= tol * norm(projWiki,np.inf)

# Testing post-processing

def test_interpolation_fits_actual_projections():
	"""
		The interpolator must fit the actual data
	"""
	res.interpolate_projection()

	T = res.projections_interpolator
	t = p.get_time_range()
	phi = p.get_alpha_range()
	Tarray = T(phi,t)

	np.testing.assert_array_almost_equal(Tarray,res.projections)

def test_fanbeam_operator_has_compact_support():
	"""
		Because we image a finite object, :math:`T(\phi/2,t)` and
		:math:`T(-\phi/2,t)` are always 0.
	"""
	res.interpolate_projection()
	T = res.projections_interpolator

	t = np.linspace(-1000,1000)

	np.testing.assert_almost_equal(T(np.pi/2, t).sum(), 0.0)
	np.testing.assert_almost_equal(T(-np.pi/2, t).sum(), 0.0)

## Testing DCC
def test_alpha_t_0_0():
	"""
		When :math:`v = 0` and :math:`x = 0`, one has :math:`\alpha(T/2) = \pi/2`
		and :math:`\alpha(-T/2) = -\pi/2`.
		Moreover, the function is not supposed to raise a warning about
		division by 0.
	"""
	alpha_t_0_0 = np.vectorize(lambda t:res.alpha(t,0,0))

	# extreme values
	np.testing.assert_almost_equal(alpha_t_0_0(p.T/2), np.pi/2)
	np.testing.assert_almost_equal(alpha_t_0_0(-p.T/2), -np.pi/2)

	# elsewhere, compare with theoretical values
	t = np.linspace(-p.T/+.01, p.T/2-.01, 1000)
	alpha_theo = np.arctan( (p.R0*np.sin(p.omega*t/360 *2*np.pi)) \
		                  / (p.R0*np.cos(p.omega*t/360 *2*np.pi)-p.y0) )
	np.testing.assert_array_almost_equal(alpha_t_0_0(t), alpha_theo)







