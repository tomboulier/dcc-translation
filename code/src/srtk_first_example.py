# SimpleRTK first example
from __future__ import print_function
import SimpleRTK as srtk
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if len(sys.argv) < 2:
    print("Usage: RTKFirstReconstruction <output>")
    sys.exit(1)

# Defines the RTK geometry object
geometry = srtk.ThreeDCircularProjectionGeometry()
numberOfProjections = 360
firstAngle = 0
angularArc = 360
sid = 600  # source to isocenter distance in mm
sdd = 1200  # source to detector distance in mm
isox = 0  # X coordinate on the projection image of isocenter
isoy = 0  # Y coordinate on the projection image of isocenter
for x in range(0, numberOfProjections):
    angle = firstAngle + x * angularArc / numberOfProjections
    geometry.AddProjection(sid, sdd, angle, isox, isoy)

constantImageSource = srtk.ConstantImageSource()
origin = [ -127.5, -127.5, 0. ]
sizeOutput = [ 256, 256,  numberOfProjections ]
spacing = [ 1.0, 1.0, 1.0 ]
constantImageSource.SetOrigin( origin )
constantImageSource.SetSpacing( spacing )
constantImageSource.SetSize( sizeOutput )
constantImageSource.SetConstant(0.0)
source = constantImageSource.Execute()
 
rei = srtk.RayEllipsoidIntersectionImageFilter()
semiprincipalaxis = [ 50, 50, 50]
center = [ 0, 0, 0]
# Set GrayScale value, axes, center...
rei.SetDensity(20)
rei.SetAngle(0)
rei.SetCenter(center)
rei.SetAxis(semiprincipalaxis)
rei.SetGeometry( geometry )
reiImage = rei.Execute(source)
 
# Create reconstructed image
constantImageSource2 = srtk.ConstantImageSource()
origin = [ -63.5, -63.5, -63.5 ]
sizeOutput = [ 128, 128, 128 ]
constantImageSource2.SetOrigin( origin )
constantImageSource2.SetSpacing( spacing )
constantImageSource2.SetSize( sizeOutput )
constantImageSource2.SetConstant(0.0)
source2 = constantImageSource2.Execute()
 
print("Performing reconstruction")
feldkamp = srtk.FDKConeBeamReconstructionFilter()
feldkamp.SetGeometry( geometry );
feldkamp.SetTruncationCorrection(0.0);
feldkamp.SetHannCutFrequency(0.0);
image = feldkamp.Execute(source2,reiImage) 
 
print("Masking field-of-view")
fov = srtk.FieldOfViewImageFilter()
fov.SetGeometry(geometry)
fov.SetProjectionsStack(reiImage)
image = fov.Execute(image)

plt.imshow(srtk.GetArrayFromImage(image[:,64,:]), cmap = cm.Greys_r)
plt.show()
 
writer = srtk.ImageFileWriter()
writer.SetFileName ( sys.argv[1] )
writer.Execute ( image );