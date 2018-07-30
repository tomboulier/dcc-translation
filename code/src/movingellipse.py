#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
    This script simulates a moving ellipse, illuminated by a source
    following an arc of circle. The illumination is fanbeam.
"""

from simulator import Simulator
from parameters import Parameters
from dcc import DataConsistencyConditions
from dcc import PolynomProjector
import math
import numpy as np
import ConfigParser
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
matplotlib.rcParams['text.latex.unicode'] = True


class Optimizer(object):
    """
        Optimization process in order to recover velocity
        occurs in this class
    """

    def __init__(self, polynom_projector):
        """
            In order to avoid inverse crime, we have to be
            very careful here with the parameters the DCCs
            can access.
        """
        self.polyproj = polynom_projector

    def minimize_rmse(self, n, x):
        """
            Minimization of RMSE of the fitting
        """
        from scipy.optimize import minimize

        residual_callable = lambda v: self.polyproj.residual_polyfit(n, v, x)
        # residual_callable = lambda v: self.polyproj.compute_RMSE(n,v,x)

        min_v = minimize(residual_callable, 0, method='Powell')

        return min_v.x


if __name__ == '__main__':
    p = Parameters('example.ini')
    res = Simulator(p).run()
    res.plotSinogram()

    # Plot DCC
    DCC = DataConsistencyConditions(res)

    polyproj = PolynomProjector(DCC)
    x = DCC.get_virtual_points_vector(p.v, p.v2)[10:40]
    polyproj.plot_fitting(1, p.v, p.v2, x)
