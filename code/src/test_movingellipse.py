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










