from scanner import Source
from scanner import Detector
from target import *
from sinogram import Results
import numpy as np

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
        if params.targetType == "ellipse":
            if params.ellipseModel == "RTK":
                self.target = MovingEllipse(params)
            elif params.ellipseModel == "Analytical":
                self.target = AnalyticalEllipse(params)
            else:
                raise ValueError("Type of Ellipse can only be 'RTK' or 'Analytical'")
        else:
            raise NotImplementedError("Simulation only works for ellipse target")

    def run(self):
        """
            Computes all the projections, at any angle and any time
        """
        # get general parameters
        T = self.params.T
        imageSize = self.params.imageSize
        Nt = self.params.Nt

        # the array giving projections
        projarray = np.zeros((Nt, imageSize))

        for nt, t in enumerate(self.params.get_time_range()):
            projarray[nt, :] = self.target.compute_projection(t,
                                                              self.source,
                                                              self.detector)

        # store results
        results = Results(self.params, self.source, self.detector)
        results.projections = projarray
        return results
