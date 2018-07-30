from main import *
from numpy.linalg import norm
import pytest
import SimpleRTK as srtk


def wiki_example_RTK(configFile):
    """
        Simulates a cone-beam illumination of an ellipse.

        Taken from RTK's wiki, see the following page:
        http://wiki.openrtk.org/index.php/SimpleRTK#Testing_SimpleRTK
    """
    # Read the parameters from file
    cfg = ConfigParser.ConfigParser()
    cfg.read(configFile)
    numberOfProjections = cfg.getint('Parameters', 'Nt')
    omega = cfg.getfloat('Parameters', 'omega')
    T = cfg.getfloat('Parameters', 'T')
    firstAngle = 90 - omega * T / 2
    angularArc = omega * T
    sid = cfg.getfloat('Parameters', 'R0')  # source to isocenter distance in mm
    sdd = cfg.getfloat('Parameters', 'sdd')  # source to detector distance in mm
    isox = 0  # X coordinate on the projection image of isocenter
    isoy = 0  # Y coordinate on the projection image of isocenter
    imageSize = cfg.getint('Parameters', 'imageSize')

    # get ellipse parameters
    ellipseDensity = cfg.getfloat('Ellipse', 'density')
    ellipseAngle = cfg.getfloat('Ellipse', 'angle')
    ellipseCenterX = cfg.getfloat('Ellipse', 'x')
    ellipseCenterY = cfg.getfloat('Ellipse', 'y')
    ellipseSemiAxisX = cfg.getfloat('Ellipse', 'a')
    ellipseSemiAxisY = cfg.getfloat('Ellipse', 'b')

    # Defines the RTK geometry object
    geometry = srtk.ThreeDCircularProjectionGeometry()
    for x in range(0, numberOfProjections):
        angle = firstAngle + x * angularArc / (numberOfProjections - 1)
        geometry.AddProjection(sid, sdd, angle, isox, isoy)

    # initialization of 'blank' image
    constantImageSource = srtk.ConstantImageSource()
    sizeOutput = [imageSize, imageSize, numberOfProjections]
    spacing = [1.0, 1.0, 1.0]
    origin = (np.array(sizeOutput) - 1) * np.array(spacing) * -.5
    constantImageSource.SetOrigin(origin)
    constantImageSource.SetSpacing(spacing)
    constantImageSource.SetSize(sizeOutput)
    constantImageSource.SetConstant(0.0)
    source = constantImageSource.Execute()

    rei = srtk.RayEllipsoidIntersectionImageFilter()
    semiprincipalaxis = [ellipseSemiAxisX,
                         ellipseSemiAxisY,
                         ellipseSemiAxisY]
    center = [ellipseCenterX, ellipseCenterY, 0]
    # Set GrayScale value, axes, center...
    rei.SetDensity(ellipseDensity)
    rei.SetAngle(ellipseAngle)
    rei.SetCenter(center)
    rei.SetAxis(semiprincipalaxis)
    rei.SetGeometry(geometry)
    reiImage = rei.Execute(source)

    return srtk.GetArrayFromImage(reiImage[:, imageSize / 2, :])


configFile = 'test.ini'
tol = .01

# simulation by our algorithm
p = Parameters(configFile)
s = Simulator(p)
res = s.run()
DCC = DataConsistencyConditions(res)

# simulation with RTK's example in wiki
projWiki = wiki_example_RTK(configFile)

# difference between the two projections
diff = res.projections - projWiki


## General tests
def test_get_angle():
    """
        The angles start at 90 since we compute the angle with the
        y-axis
    """
    assert p.get_angle(0) == 90


def test_y0():
    """
        :math:`y0` is supposed to be stricly positive
    """
    assert p.y0 > 0


## Testing projection operator

def test_projections_have_proper_dimensions():
    assert res.projections.shape[0] == p.Nt
    assert res.projections.shape[1] == p.imageSize


def test_projections_fits_example_in_wiki_RTK():
    """
        Compare the two projections, in L2 norm
    """
    assert norm(diff) <= tol * norm(projWiki)
    assert norm(diff, 1) <= tol * norm(projWiki, 1)
    assert norm(diff, 2) <= tol * norm(projWiki, 2)
    assert norm(diff, np.inf) <= tol * norm(projWiki, np.inf)


## Testing post-processing

configFile = 'example.ini'
tol = .01

# simulation by our algorithm
p = Parameters(configFile)
s = Simulator(p)
res = s.run()
DCC = DataConsistencyConditions(res)


def test_interpolation_fits_actual_projections():
    """
        The interpolator must fit the actual data
    """
    res.interpolate_projection()

    T = res.projections_interpolator
    t = p.get_time_range()
    phi = p.get_alpha_range()
    Tarray = T(phi, t)

    np.testing.assert_array_almost_equal(Tarray, res.projections)


def test_fanbeam_operator_has_compact_support():
    """
        Because we image a finite object, :math:`T(\phi/2,t)` and
        :math:`T(-\phi/2,t)` are always 0.
    """
    res.interpolate_projection()
    T = res.projections_interpolator

    t = np.linspace(-1000, 1000)

    np.testing.assert_almost_equal(T(np.pi / 2, t).sum(), 0.0)
    np.testing.assert_almost_equal(T(-np.pi / 2, t).sum(), 0.0)


## Testing DCC
def test_lambda_t_0_0():
    """
        When :math:`v = 0` and :math:`x = 0`, one has :math:`\lambda(T/2) = \pi/2`
        and :math:`\lambda(-T/2) = -\pi/2`.
        Moreover, the function is not supposed to raise a warning about
        division by 0.
    """
    lambda_t_0_0 = np.vectorize(lambda t: DCC.lambda_v(t, 0, 0, 0))

    # extreme values
    np.testing.assert_almost_equal(lambda_t_0_0(p.T / 2), np.pi / 2)
    np.testing.assert_almost_equal(lambda_t_0_0(-p.T / 2), -np.pi / 2)

    # elsewhere, compare with theoretical values
    t = np.linspace(-p.T / +.01, p.T / 2 - .01, 1000)
    lambda_theo = np.arctan((p.R0 * np.sin(p.omega * t / 360 * 2 * np.pi)) \
                            / (p.R0 * np.cos(p.omega * t / 360 * 2 * np.pi) - p.y0))
    np.testing.assert_array_almost_equal(lambda_t_0_0(t), lambda_theo)


def test_beta_v2_0():
    """
        When :math:`v_2 = 0`, then :math:`\beta = 0`.
        Here, we test that beta(v1,0) is equal to 0, with a random
        value (strictly positive) for v1.
    """
    from numpy.random import randn
    v1 = np.abs(randn())

    assert DCC.beta(v1, 0) == 0.


def test_beta_raises_exception():
    """
        Assert that method beta in DCC class raises the proper
        exception, ie when denominator is equal to 0.
    """
    v1 = -2 * p.R0 * np.sin(2 * np.pi / 360 * p.omega * p.T / 2) / p.T

    with pytest.raises(ZeroDivisionError) as e_info:
        DCC.beta(v1, p.v2)


def test_virtual_source_position_begin_end():
    """
        By definition, the y-coordinates of the source is supposed
        to be the same when :math:`t=-T/2` and :math:`t=T/2`.
    """
    # y-coordinate of s_v(T/2)
    s_v_2_T_over_2 = DCC.get_virtual_source_position(p.T / 2, p.v, p.v2)[1]

    # y-coordinate of s_v(T/2)
    s_v_2_minus_T_over_2 = DCC.get_virtual_source_position(-p.T / 2, p.v, p.v2)[1]

    np.testing.assert_almost_equal(s_v_2_T_over_2, s_v_2_minus_T_over_2)


def test_simplified_formula_for_y0_prime():
    """
        Test that :math:`s^{(v)}(T/2) = s^{(v)}(-T/2) = R_0 \cos(\omega T/2 + \beta)`
    """
    y0_prime = p.R0 * np.cos(p.omega / 360 * 2 * np.pi * p.T / 2 + DCC.beta(p.v, p.v2))
    # y-coordinate of s_v(T/2)
    s_v_2_T_over_2 = DCC.get_virtual_source_position(p.T / 2, p.v, p.v2)[1]

    # y-coordinate of s_v(T/2)
    s_v_2_minus_T_over_2 = DCC.get_virtual_source_position(-p.T / 2, p.v, p.v2)[1]

    np.testing.assert_almost_equal(y0_prime, s_v_2_minus_T_over_2)
    np.testing.assert_almost_equal(s_v_2_T_over_2, y0_prime)


def test_F_0_0():
    """
        As mentioned in equation (17) of the abstract, when :math:`v=0`
        the function :math:`F` has a simpler formulation

        ..math::
        F^{(0,0)}(t,x) = \frac{x + R_0 \sin(\omega t)}{R0 \cos(\omega t)}
    """
    R0 = p.R0
    omega = p.omega / 360 * 2 * np.pi
    y0 = R0 * np.cos(omega * p.T / 2)

    from numpy.random import randn
    # virtual point, taken randomly
    x = randn()
    # time vector, taken randomly
    t = randn()

    F_theo = (x + R0 * np.sin(omega * t)) / (R0 * np.cos(omega * t) - y0)
    np.testing.assert_almost_equal(DCC.F(t, 0, 0, x), F_theo)


def test_get_virtual_source_position_v_0():
    """
        Test that when :math:`v_1 = v_2 = 0`, one has
        :math:`s^{(v)}(t) = R_0(-\sin(\omega t), \cos(\omega t))`
    """
    from numpy.random import randn
    t = randn()

    s_0_t = DCC.get_virtual_source_position(t, 0, 0)
    omega = p.omega / 360 * 2 * np.pi
    s_theo = p.R0 * np.array((-np.sin(omega * t), np.cos(omega * t)))

    np.testing.assert_almost_equal(s_0_t, s_theo)


def test_get_virtual_source_position_v2_0():
    """
        Test that when :math:`v_2 = 0`, one has
        :math:`s^{(v)}(t) = (-R_0 \sin(\omega t) - (t+T/2) v_1, \
                              R_0 \cos(\omega t))`
    """
    from numpy.random import randn
    t = randn()
    v = randn()
    R0 = p.R0
    omega = p.omega / 360 * 2 * np.pi
    T = p.T

    s_0_t = DCC.get_virtual_source_position(t, v, 0)
    omega = p.omega / 360 * 2 * np.pi

    # we will use here the same code as in previous versions
    Mvt = np.array(((t + p.T / 2) * v, 0))
    s_theo = s.source.get_position(t) - Mvt

    s_theo_analytical = np.array((-R0 * np.sin(omega * t) - \
                                  (t + T / 2) * v, R0 * np.cos(omega * t)))

    np.testing.assert_almost_equal(s_0_t, s_theo)
    np.testing.assert_almost_equal(s_0_t, s_theo_analytical)


def test_get_virtual_points_vector():
    """
        Here, we test that when :math:`v_1 = v_2 = 0`, the virtual points
        :math:`x'` are in fact points between :math:`-R_0\sin(\omega T/2)`
        and :math:`-R_0\sin(\omega T/2)`
    """
    x_prime = DCC.get_virtual_points_vector(0, 0)

    omega = p.omega / 360 * 2 * np.pi
    x_prime_theo = p.R0 * np.linspace(-np.sin(omega * p.T / 2), np.sin(omega * p.T / 2))

    np.testing.assert_almost_equal(x_prime, x_prime_theo)


def test_jacobian_with_old_formulas():
    """
        In some archives, we had the formulas for the jacobian
        when :math:`v_2=0`, so we will use it here.
    """

    # for time, degree of polynom, position: we take random values
    from numpy.random import randint
    n = randint(5)
    from numpy.random import randn
    t = np.abs(randn())
    x = randn()
    v1 = np.abs(randn())
    # important, because the test relies on the fact that v2 = 0
    v2 = 0

    # usual parameters
    R0 = p.R0
    T = p.T
    omega = p.omega / 360 * 2 * np.pi
    y0 = p.y0

    # compute weighting function with DCC class
    jacobian_new = DCC.jacobian(t, v1, v2, x)

    # compute old formula (with previous version of code)
    v = v1  # to be compatible with old notations
    # first, "jacobian"
    num = R0 ** 2 * omega - v * y0 + R0 * np.cos(omega * t) * (v - omega * y0) \
          + R0 * omega * np.sin(omega * t) * (x + (t + T / 2) * v)
    denom = (R0 * np.cos(omega * t) - y0) ** 2
    jacobian_old = num / denom

    np.testing.assert_almost_equal(jacobian_new, jacobian_old)


def test_integrand_with_old_formulas():
    """
        In some archives, we had the formulas for the weighting function
        and the jacobian when :math:`v_2=0`, so we will use it here.
    """

    # for time, degree of polynom, position: we take random values
    from numpy.random import randint
    n = randint(5)
    from numpy.random import randn
    t = np.abs(randn())
    x = randn()
    v1 = np.abs(randn())
    # important, because the test relies on the fact that v2 = 0
    v2 = 0

    # usual parameters
    R0 = p.R0
    T = p.T
    omega = p.omega / 360 * 2 * np.pi
    y0 = R0 * np.cos(omega * T / 2)

    # compute weighting function with DCC class
    weight_new = DCC.W(n, t, v1, v2, x)

    # compute old formula (with previous version of code)
    v = v1  # to be compatible with old notations
    # first, "jacobian"
    num = R0 ** 2 * omega - v * y0 + R0 * np.cos(omega * t) * (v - omega * y0) \
          + R0 * omega * np.sin(omega * t) * (x + (t + T / 2) * v)
    denom = (R0 * np.cos(omega * t) - y0) ** 2
    J_x_t_v = num / denom

    # then, function W_n(t,x)
    s_v_t = DCC.get_virtual_source_position(t, v1, v2)
    from scipy.spatial import distance
    D_x_t = distance.euclidean(s_v_t, (x, y0))
    num = (x + R0 * np.sin(omega * t) + (t + T / 2) * v) ** n
    denom = D_x_t * (R0 * np.cos(omega * t) - y0) ** (n - 1)
    W_n_t_x = num / denom

    weight_old = J_x_t_v * W_n_t_x

    np.testing.assert_almost_equal(weight_old, weight_new)


def test_two_different_formulas_for_weight_function():
    """
        In Theorem 1, the weight function is defined as

        ..math::
        W(t) = \tan^n(\lambda) \cos(\lambda) \partial_t \lambda

        We can also show that

        ..math::
        W(t) = \frac{F(t)^n}{\sqrt{1+F^(t)}} F'(t),

        where :math:`F(t) := \arctan(\lambda)`. This is what
        we will test here.
    """
    # for time, degree of polynom, position, velocity,
    # we take random values
    from numpy.random import randint
    n = randint(5)
    from numpy.random import randn
    t = np.abs(randn())
    x = randn()
    v1 = np.abs(randn())
    v2 = np.abs(randn())

    # usual parameters
    R0 = p.R0
    T = p.T
    omega = p.omega / 360 * 2 * np.pi
    y0 = R0 * np.cos(omega * T / 2)

    # compute weighting function with DCC class
    weight_DCC = DCC.W(n, t, v1, v2, x)

    # compute formula
    F_t_x = DCC.F(t, v1, v2, x)
    diff_t_F = DCC.jacobian(t, v1, v2, x)  # computes F'(t)
    weight_alternative = F_t_x ** n / np.sqrt(1 + F_t_x ** 2) * diff_t_F

    np.testing.assert_almost_equal(weight_DCC, weight_alternative)
