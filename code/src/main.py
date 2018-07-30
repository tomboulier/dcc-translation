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
from inverseproblem import Optimizer
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
matplotlib.rcParams['text.latex.unicode'] = True


if __name__ == '__main__':
    p = Parameters('example.ini')
    res = Simulator(p).run()
    res.plotSinogram()

    # Plot DCC
    DCC = DataConsistencyConditions(res)

    polyproj = PolynomProjector(DCC)
    n=1
    x = DCC.get_virtual_points_vector(p.v, p.v2)[10:40]
    polyproj.plot_fitting(n, p.v, p.v2, x)

    # # recover velocity
    # optimizer = Optimizer(polyproj)
    # v_optim = optimizer.minimize_rmse(n, x)
    #
    # print "Result of optimization is : " + str(v_optim)
