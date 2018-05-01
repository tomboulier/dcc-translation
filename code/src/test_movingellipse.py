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
DCC = DataConsistencyConditions(res)

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

## Testing post-processing

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
def test_lambda_t_0_0():
	"""
		When :math:`v = 0` and :math:`x = 0`, one has :math:`\lambda(T/2) = \pi/2`
		and :math:`\lambda(-T/2) = -\pi/2`.
		Moreover, the function is not supposed to raise a warning about
		division by 0.
	"""
	lambda_t_0_0 = np.vectorize(lambda t:DCC.lambda_v(t,0,0,0))

	# extreme values
	np.testing.assert_almost_equal(lambda_t_0_0(p.T/2), np.pi/2)
	np.testing.assert_almost_equal(lambda_t_0_0(-p.T/2), -np.pi/2)

	# elsewhere, compare with theoretical values
	t = np.linspace(-p.T/+.01, p.T/2-.01, 1000)
	lambda_theo = np.arctan( (p.R0*np.sin(p.omega*t/360 *2*np.pi)) \
		                  / (p.R0*np.cos(p.omega*t/360 *2*np.pi)-p.y0) )
	np.testing.assert_array_almost_equal(lambda_t_0_0(t), lambda_theo)

def test_beta_v2_0():
	"""
		When :math:`v_2 = 0`, then :math:`\beta = 0`.
		Here, we test that beta(v1,0) is equal to 0, with a random
		value (strictly positive) for v1.
	"""
	from numpy.random import randn
	v1 = np.abs(randn())

	assert DCC.beta(v1,0) == 0.

def test_beta_raises_exception():
	"""
		Assert that method beta in DCC class raises the proper
		exception, ie when denominator is equal to 0.
	"""
	v1 = -2*p.R0*np.sin(2*np.pi/360 * p.omega*p.T/2)/p.T

	with pytest.raises(ZeroDivisionError) as e_info:
		DCC.beta(v1,p.v2)

def test_virtual_source_position_begin_end():
	"""
		By definition, the y-coordinates of the source is supposed
		to be the same when :math:`t=-T/2` and :math:`t=T/2`.
	"""
	# y-coordinate of s_v(T/2)
	s_v_2_T_over_2 = DCC.get_virtual_source_position(p.T/2, p.v, p.v2)[1]

	# y-coordinate of s_v(T/2)
	s_v_2_minus_T_over_2 = DCC.get_virtual_source_position(-p.T/2, p.v, p.v2)[1]
	

	np.testing.assert_almost_equal(s_v_2_T_over_2, s_v_2_minus_T_over_2)

def test_simplified_formula_for_y0_prime():
	"""
		Test that :math:`s^{(v)}(T/2) = s^{(v)}(-T/2) = R_0 \cos(\omega T/2 + \beta)`
	"""
	y0_prime = p.R0 * np.cos(p.omega/360*2*np.pi * p.T/2 + DCC.beta(p.v,p.v2))
	# y-coordinate of s_v(T/2)
	s_v_2_T_over_2 = DCC.get_virtual_source_position(p.T/2, p.v, p.v2)[1]

	# y-coordinate of s_v(T/2)
	s_v_2_minus_T_over_2 = DCC.get_virtual_source_position(-p.T/2, p.v, p.v2)[1]

	np.testing.assert_almost_equal(y0_prime, s_v_2_minus_T_over_2)
	np.testing.assert_almost_equal(s_v_2_T_over_2, y0_prime)

def test_F_0_0():
	"""
		As mentioned in equation (17) of the abstract, when :math:`v=0`
		the function :math:`F` has a simpler formulation

		..math::
		F^{(0,0)}(t,x) = \frac{x + R_0 \sin(\omega t)}{R0 \cos(\omega t)}
	"""
	R0 = p.R0
	omega = p.omega/360 *2*np.pi
	y0 = R0*np.cos(omega*p.T/2)
	
	from numpy.random import randn
	# virtual point, taken randomly
	x = randn()
	# time vector, taken randomly
	t = randn()
	
	F_theo = (x + R0 * np.sin(omega * t)) / (R0 * np.cos(omega*t) - y0)
	np.testing.assert_almost_equal(DCC.F(t, 0, 0, x), F_theo)

def test_get_virtual_source_position():
	"""
		Test that when :math:`v_1 = v_2 = 0`, one has
		:math:`s^{(v)}(t) = R_0(-\sin(\omega t), \cos(\omega t))`
	"""
	from numpy.random import randn
	t = randn()

	s_0_t = DCC.get_virtual_source_position(t, 0, 0)
	omega = p.omega/360*2*np.pi
	s_theo = p.R0 * np.array( (-np.sin(omega*t), np.cos(omega*t)) )

	np.testing.assert_almost_equal(s_0_t, s_theo)

def test_get_virtual_points_vector():
	"""
		Here, we test that when :math:`v_1 = v_2 = 0`, the virtual points
		:math:`x'` are in fact points between :math:`-R_0\sin(\omega T/2)`
		and :math:`-R_0\sin(\omega T/2)`
	"""
	x_prime = DCC.get_virtual_points_vector(0,0)

	omega = p.omega/360*2*np.pi
	x_prime_theo = p.R0 * np.linspace(-np.sin(omega*p.T/2), np.sin(omega*p.T/2) )

	np.testing.assert_almost_equal(x_prime, x_prime_theo)



