import numpy as np
import ConfigParser

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
        self.targetType = cfg.get('Parameters', 'targetType')

        if self.targetType == "ellipse":
            # parameters describing the ellipse
            self.ellipseDensity = cfg.getfloat('Ellipse', 'density')
            self.ellipseAngle = cfg.getfloat('Ellipse', 'angle')
            self.ellipseCenterX = cfg.getfloat('Ellipse', 'x')
            self.ellipseCenterY = cfg.getfloat('Ellipse', 'y')
            self.ellipseSemiAxisX = cfg.getfloat('Ellipse', 'a')
            self.ellipseSemiAxisY = cfg.getfloat('Ellipse', 'b')
            self.ellipseModel = cfg.get('Ellipse', 'model')

        # These parameters are derived from previous ones
        self.y0 = self.R0 * np.cos(self.omega * self.T / 2 / 360 * 2 * np.pi)

    def get_angle(self, t, unit='degrees'):
        return self.omega * t + 90

    def get_max_angle(self):
        return self.get_angle(self.T / 2)

    def get_min_angle(self):
        return self.get_angle(-self.T / 2)

    def get_time_range(self):
        nt = np.arange(self.Nt)
        return -self.T / 2 + self.T * nt / (self.Nt - 1)

    def get_alpha_range(self):
        x = np.arange(-self.imageSize / 2, self.imageSize / 2)
        # return np.arctan( x / self.sdd) * 360 / (2*np.pi)
        return np.arctan(x / self.sdd)
