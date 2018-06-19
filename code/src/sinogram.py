import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate
plt.rc('text', usetex=True)
matplotlib.rcParams['text.latex.unicode'] = True

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

    def plotSinogram(self, xunits='mm'):
        # define the limits of the axis
        imageSize = self.params.imageSize
        T = self.params.T
        max_angle = self.params.get_max_angle()
        min_angle = self.params.get_min_angle()
        phimax = np.arctan(.5 * imageSize / self.params.sdd) * 360 / (2 * np.pi)

        # plot the image
        plt.figure()

        if xunits == 'mm':
            # the units here represent a distance (on the detector)
            plt.xlabel('Distance from detector center (in mm)', labelpad=20)
            extent = [-imageSize / 2, imageSize / 2, max_angle, min_angle]
            aspect = imageSize / (max_angle - min_angle)

        elif xunits == 'degrees':
            # the units here represent an angle ('phi' in T(x,phi))
            plt.xlabel('Beam direction (in degrees)', labelpad=20)
            extent = [-phimax, phimax, max_angle, min_angle]
            aspect = 2 * phimax / (max_angle - min_angle)

        plt.imshow(self.projections, cmap=cm.Greys_r, extent=extent,
                   aspect=aspect / 2)
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
