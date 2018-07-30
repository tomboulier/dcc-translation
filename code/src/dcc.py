import numpy as np
import math
from scipy import integrate
import matplotlib.pyplot as plt
import matplotlib

class DataConsistencyConditions(object):
    """
        The computations of DCCs happen here
    """

    def __init__(self, results):
        self.results = results
        self.params = results.params
        self.source = results.source

        self.DCC_function = None

    def get_virtual_source_position(self, t, v1, v2):
        """
            Computes :math:`s_v(t)` given by

            ..math::
            s_{\mathbf{v}}(t) = \mathcal{R}_{-\beta} \left( s(t) - M_v(t) \right)
        """
        # translation vector
        Mvt = (t + self.params.T / 2) * np.array((v1, v2))

        # rotation matrix
        beta = self.beta(v1, v2)
        c, s = np.cos(-beta), np.sin(-beta)
        R_beta = np.array(((c, -s), (s, c)))

        return R_beta.dot(self.source.get_position(t) - Mvt)

    def get_virtual_points_vector(self, v1, v2):
        """
            These will be the x-coordinates of points where DCC function
            will be computed.
        """
        # xmax = self.R0 * np.sin(self.omega*self.T/2 /360*2*np.pi) * .75
        # return np.linspace(-xmax, xmax)

        # take notations of Definition 3 for x', i.e. x' is between x'_l and x'_r
        # whith x'_l = s_v(T/2) and x'_r = s_v(-T/2)
        T = self.params.T
        # v1 = self.params.v
        # v2 = self.params.v2

        xl = self.get_virtual_source_position(T / 2, v1, v2)[0]
        xr = self.get_virtual_source_position(-T / 2, v1, v2)[0]

        return np.linspace(xl, xr)

    def beta(self, v1, v2):
        """
            The rotation angle :math:`\beta` defined in formula (7)
            of the abstract

            .. math::
            \beta = \arctan \left( \frac{T v_2} \
            {2 R_0 \sin(\omega T/2) + T v_1} \right)
        """
        R0 = self.params.R0
        omega = self.params.omega / 360 * 2 * np.pi
        T = self.params.T

        num = T * v2
        denom = 2 * R0 * np.sin(omega * T / 2) + T * v1

        if denom == 0.:
            raise ZeroDivisionError("denominator is equal to zero in function beta")
        else:
            return np.arctan(num / denom)

    def lambda_v(self, t, v1, v2, x):
        """
            Function :math:`\lambda^{(v)}` in Theorem 1, i.e.

            ..math::
            \lambda^{(\bv)}(t,x') = \arctan \left( F^{(\bv)}(t,x') \right),

            where :math:`F^{(\bv)}(t,x')` is defined as the fraction
            :math:`A/B` with

            ..math::
                A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) + \
                \sin \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right)

            and
            ..math::
                B = \cos \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right) - \
                \sin \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) - y'_0

        """

        return np.arctan(self.F(t, v1, v2, x))

    def F(self, t, v1, v2, x):
        """
            Function :math:`F^{(\bv)}(t,x')` which is defined in
            Theorem 1 as the fraction :math:`A/B` with

            ..math::
                A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) + \
                \sin \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right)

            and
            ..math::
                B = \cos \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right) - \
                \sin \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) - y'_0

        """

        A = self.F_numerator(t, v1, v2, x)
        B = self.F_denominator(t, v1, v2, x)

        if B == 0.:
            return np.sign(A) * np.inf
        else:
            return A / B

    def F_numerator(self, t, v1, v2, x):
        """
            The numerator of the function :math:`F^{(\bv)}(t,x')`
            defined in Theorem, where it is denoted :math:`A`, with

            ..math::
                A = x' + \cos \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) + \
                \sin \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right)
        """
        R0 = self.params.R0
        omega = self.params.omega / 360 * 2 * np.pi
        T = self.params.T
        beta = self.beta(v1, v2)

        return x + np.cos(beta) * (R0 * np.sin(omega * t) + (t + T / 2) * v1) + \
               np.sin(beta) * (R0 * np.cos(omega * t) - (t + T / 2) * v2)

    def F_denominator(self, t, v1, v2, x):
        """
            The denominator of the function :math:`F^{(\bv)}(t,x')`
            defined in Theorem, where it is denoted :math:`B`, with

            ..math::
                B = \cos \beta \left( R_0 \cos(\omega t) - \
                \left( t + T/2 \right)v_2 \right) - \
                \sin \beta \left( R_0 \sin(\omega t) + \
                \left( t + T/2 \right)v_1 \right) - y'_0
        """
        R0 = self.params.R0
        omega = self.params.omega / 360 * 2 * np.pi
        T = self.params.T
        beta = self.beta(v1, v2)
        y0_prime = R0 * np.cos(omega * T / 2 + beta)

        return np.cos(beta) * (R0 * np.cos(omega * t) - (t + T / 2) * v2) - \
               np.sin(beta) * (R0 * np.sin(omega * t) + (t + T / 2) * v1) - y0_prime

    def diff_t_F_numerator(self, t, v1, v2, x):
        """
            Derivative of the function :math:`A` defined in equation (13),
            with respect to time.
        """
        R0 = self.params.R0
        omega = self.params.omega / 360 * 2 * np.pi
        T = self.params.T
        beta = self.beta(v1, v2)

        return np.cos(beta) * (R0 * omega * np.cos(omega * t) + v1) + \
               np.sin(beta) * (- R0 * omega * np.sin(omega * t) - v2)

    def diff_t_F_denominator(self, t, v1, v2, x):
        """
            Derivative of the function :math:`B` defined in equation (14),
            with respect to time.
        """
        R0 = self.params.R0
        omega = self.params.omega / 360 * 2 * np.pi
        T = self.params.T
        beta = self.beta(v1, v2)

        return np.cos(beta) * (- R0 * omega * np.sin(omega * t) - v2) - \
               np.sin(beta) * (R0 * omega * np.cos(omega * t) + v1)

    def jacobian(self, t, v1, v2, x):
        """
            The Jacobian used in the change of variable in Theorem 1,
            which is given by the derivative of :math:`\lambda^{(\mathbf{v})}(t,x')`
            with respect to :math:`t`. Hence, one has

            ..math::
            \partial_t \lambda^{(\mathbf{v})} = \
            \frac{\partial_t A B - A \partial_t}{B^2} \frac{1}{1 + F^2}
        """
        A = self.F_numerator(t, v1, v2, x)
        B = self.F_denominator(t, v1, v2, x)
        diff_t_A = self.diff_t_F_numerator(t, v1, v2, x)
        diff_t_B = self.diff_t_F_denominator(t, v1, v2, x)
        F = self.F(t, v1, v2, x)

        # return (diff_t_A * B - A * diff_t_B) / (B*B * (1+F*F))
        return (diff_t_A * B - A * diff_t_B) / (B * B)

    def W(self, n, t, v1, v2, x):
        """
            The weighting function in Theorem 1, i.e.

            ..math::
            W_n^{(v)}(t,x') = \tan^n \left( \lambda^{(v)}(t,x') \right) \\
            \cos \left( \lambda^{(v)}(t,x') \right) J(x,t,v),

            where :math:`J(x,t,v)` is given by the function 'jacobian'
        """
        lambda_v = self.lambda_v(t, v1, v2, x)
        jacobian = self.jacobian(t, v1, v2, x)

        return (np.tan(lambda_v)) ** n * np.cos(lambda_v) * jacobian

    def integrand(self, n, t, v1, v2, x):
        """
            The function to be integrated in DCC, i.e.

            ..math::
            \mathcal{F}(t,\lambda^{(v)}(t,x)) W_n^{(v)}(t,x)

            We put it in a function since we have to be careful with Nans
        """

        T = self.params.T
        omega = self.params.omega / 360 * 2 * np.pi
        beta = self.beta(v1, v2)

        lambda_v = lambda t, x: self.lambda_v(t, v1, v2, x)
        fb_proj = lambda t, x: self.results.projections_interpolator(lambda_v(t, x) + beta - omega * t, t)
        weight = lambda t, x, n: self.W(n, t, v1, v2, x)

        y = fb_proj(t, x) * weight(t, x, n)

        if fb_proj(t, x) < 1e-10:
            return 0.

        if math.isnan(y):
            return 0.

        return y

    def B(self, n, v1, v2, x):
        """
            Compute the function B_n(x) in Theorem 1, which is supposed
            to be a polynom of order at most n, where n is the order of DCC.
        """
        T = self.params.T
        Nt = self.params.Nt

        integrand = lambda n, t, x: self.integrand(n, t, v1, v2, x)

        t = np.linspace(-T / 2, T / 2, Nt)
        y = np.array([integrand(n, time, x) for time in t])

        return integrate.simps(y, dx=t[1] - t[0])

    def compute_DCC_function(self, n, v1, v2):
        """
            Transform B as a lambda function
        """
        if self.results.projections_interpolator is None:
            self.results.interpolate_projection()

        self.DCC_function = lambda n, x: self.B(n, v1, v2, x)

    def compute_vectorized_DCC_function(self, n, v1, v2, x):
        """
            Compute a vector giving all values of B(x) for each
            point in x
        """
        self.compute_DCC_function(n, v1, v2)
        Bn = np.vectorize(lambda x: self.DCC_function(n, x))

        return Bn(x)


class PolynomProjector(object):
    """
        In this class, we fit the DCC function onto a polynomial
        function
    """

    def __init__(self, dcc):
        self.dcc = dcc

    def fit_dcc_polynom(self, n, v1, v2, x):
        """
            Interpolates the DCC function with a polynom of
            degree n
        """

        # computes the array with Bn(x) values
        y = self.dcc.compute_vectorized_DCC_function(n, v1, v2, x)

        # interpolation with polynom
        poly = np.polyfit(x, y, n)

        # the results is an array of values
        return np.poly1d(poly)(x)

    def compute_RMSE(self, n, v1, v2, x):
        """
            Compute the root mean square error of the fitting
            with fit_dcc_polynom
        """

        y = self.dcc.compute_vectorized_DCC_function(n, v1, v2, x)
        yfit = self.fit_dcc_polynom(n, v1, v2, x)

        diff = y - yfit

        return np.sqrt((diff * diff).sum()) / np.sqrt((y * y).sum())

    def plot_fitting(self, n, v1, v2, x):
        """
            Plot the result of the fitting procedure, showing
            the function Bn(x) with the polynom and the RMSE
        """
        # text for RMSE
        rmse = self.compute_RMSE(n, v1, v2, x)
        textrmse = r"$%.4f$" % (rmse)
        textstr = r"$RMSE = $" + textrmse

        y = self.dcc.compute_vectorized_DCC_function(n, v1, v2, x)
        yfit = self.fit_dcc_polynom(n, v1, v2, x)

        fig, ax = plt.subplots(1)
        plt.plot(x, y, 'ob')
        plt.plot(x, yfit, '-r')
        axes = plt.gca()
        # axes.set_ylim([y.mean()-20,y.mean()+20])
        matplotlib.rcParams.update({'font.size': 25})
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=30,
                verticalalignment='top')
        # plt.savefig('B' + str(n) + '.eps')
        plt.show()

    def residual_polyfit(self, n, v1, v2, x):
        """
            Computes the function |Bn(x)-P(x)|, where P(x) is the
            polynom obtained by fitting.
        """
        Bn = self.dcc.compute_vectorized_DCC_function(n, v1, v2, x)
        _, res, _, _, _ = np.polyfit(x, Bn, n, full=True)

        return res[0]
