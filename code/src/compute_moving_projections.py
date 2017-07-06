# This script simulate a moving object, illuminated by a source
# following an arc of circle.
# The illumination is fanbeam.

import SimpleRTK as srtk
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# parameters
R0 = 600 # radius of the circle, in mm
sdd = 2*R0 # source to detector distance, in mm
omega = 1 # speed of rotation, in degree/s
T = 90/omega # total time of experiment, in s
y0 = R0*np.cos(omega*T/2 * 2*np.pi/360) # distance from isocenter to the chord of
                                        # the arc of circle, in mm
imageSize = 256
Nt = 100 # number of discretization points for the time
v = 0 # speed of the object along the x-axis, in mm/s
projarray = np.zeros((Nt,1,imageSize)) # the array giving projections

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
	proj.SetOrigin((np.array(proj.GetSize())-1)*np.array(proj.GetSpacing())*-0.5)
	proj.SetConstant(0.0)
	source = proj.Execute()

	# compute intersection with ellipse
	rei = srtk.RayEllipsoidIntersectionImageFilter()
	rei.SetDensity(20)
	rei.SetAngle(0)
	rei.SetCenter([(t+T/2)*v, 0, 0])
	rei.SetAxis([25, 50, 50])
	rei.SetGeometry(geometry)
	reiImage = rei.Execute(source)

	# put this projection into projarray
	projarray[nt,0,:] = srtk.GetArrayFromImage(reiImage)[0,0,:]

# plot the resulting sinogram
plt.imshow(projarray[:,0,:], cmap = cm.Greys_r)
plt.show()