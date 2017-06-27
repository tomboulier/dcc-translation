#!/usr/bin/env python

# This example illustrates how to do a FBP fan-beam reconstruction using RTK.
# Since RTK only has 3D (back-projector), we create a 2D projection from the
# fan-beam projection. The pixels values of the first and the third rows are 
# actually useless but required to let the backprojection work. 
# Taken from http://wiki.openrtk.org/index.php/FanBeam
 
import SimpleRTK as srtk
import numpy as np
 
# Create a geometry
g = srtk.ThreeDCircularProjectionGeometry()
for i in range(720):
    g.AddProjection(1000,1500,i*.5)
 
# Simulate fan-beam acquisition
proj  = srtk.ConstantImageSource()
proj.SetSize([1024,1,720])
proj.SetSpacing([0.5]*3)
proj.SetOrigin((np.array(proj.GetSize())-1)*np.array(proj.GetSpacing())*-0.5)
proj = proj.Execute()
sl = srtk.SheppLoganPhantomFilter()
sl.SetGeometry(g)
proj = sl.Execute(proj)
 
writer = srtk.ImageFileWriter()
writer.SetFileName('fanbeam.tiff')
writer.Execute(proj)
 
# RTK currently has 3D (back-)projector. To overcome this, we add two rows
projarray = srtk.GetArrayFromImage(proj)
projarray = np.concatenate([projarray,projarray,projarray], 1)
projarray[:,0,:] = 1e10 # Just to illustrate that this row values are not used
projarray[:,2,:] = 1e10 # Just to illustrate that this row values are not used
proj2D = srtk.GetImageFromArray(projarray)
proj2D.SetSpacing(proj.GetSpacing())
proj2D.SetOrigin((np.array(proj2D.GetSize())-1)*np.array(proj2D.GetSpacing())*-0.5)
writer.SetFileName('conebeam.tiff')
writer.Execute(proj2D)
 
# Reconstruct using FDK
recon  = srtk.ConstantImageSource()
recon.SetSize([512,1,512])
recon.SetSpacing([0.5]*3)
recon.SetOrigin((np.array(recon.GetSize())-1)*np.array(recon.GetSpacing())*-0.5)
recon = recon.Execute()
fdk = srtk.FDKConeBeamReconstructionFilter()
fdk.SetGeometry(g)
recon = fdk.Execute(recon,proj2D)
writer.SetFileName('fdk.tiff')
writer.Execute(recon)