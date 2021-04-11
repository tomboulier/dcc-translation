from itk import RTK as rtk
import numpy as np

class Detector(object):
    """docstring for Detector"""

    def __init__(self, params):
        self.params = params

    def get_empty_image(self):
        proj = rtk.ConstantImageSource()
        proj.SetSize([self.params.imageSize, 1, 1])
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

    def get_geometry(self, t):
        geometry = rtk.ThreeDCircularProjectionGeometry()
        geometry.AddProjection(self.R0, self.sdd,
                               self.get_angle(t), 0, 0)
        return geometry

    def get_position(self, t):
        return np.array((-self.R0 * np.sin(self.omega * t / 360 * 2 * np.pi), \
                         self.R0 * np.cos(self.omega * t / 360 * 2 * np.pi)))
