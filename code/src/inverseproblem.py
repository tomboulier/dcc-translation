import warnings

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
        warnings.warn("Recovery of velocity vector (v1,v2) is implemented for v2 = 0 only")

        from scipy.optimize import minimize

        residual_callable = lambda v1: self.polyproj.residual_polyfit(n, v1, 0, x)
        # residual_callable = lambda v: self.polyproj.compute_RMSE(n,v,x)

        min_v = minimize(residual_callable, 0, method='Powell')

        return min_v.x

