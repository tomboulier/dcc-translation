from scanner import Source
from scanner import Detector
from target import MovingEllipse
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
        self.ellipse = MovingEllipse(params)

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
            projarray[nt, :] = self.ellipse.compute_projection(t,
                                                               self.source,
                                                               self.detector)

        # store results
        results = Results(self.params, self.source, self.detector)
        results.projections = projarray
        return results
