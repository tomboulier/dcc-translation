import numpy as np

class RTKEllipse(object):
    """
        Class of ellipse where projections are simulated with the module 'RTK', by Simon Rit
    """

    def __init__(self, params):
        self.params = params

    def get_density(self):
        return self.params.ellipseDensity

    def get_angle(self):
        return self.params.ellipseAngle

    def get_center(self, t):
        """
            Since it is moving, the position depends on t
        """
        T = self.params.T
        v = self.params.v
        v2 = self.params.v2

        return [self.params.ellipseCenterX - (t + T / 2) * v,
                self.params.ellipseCenterY - (t + T / 2) * v2,
                0]

    def get_axis(self):
        return [self.params.ellipseSemiAxisX,
                self.params.ellipseSemiAxisY,
                self.params.ellipseSemiAxisY]

    def compute_projection(self, t, source, detector):
        """
            Simulate fan-beam acquisition of the object with given
            source and detector, at time t
        """
        import SimpleRTK as srtk
        # create geometry of the source at time t
        geometry = source.get_geometry(t)

        # compute intersection of fan-beam with ellipse
        empty_image_detector = detector.get_empty_image()
        rei = srtk.RayEllipsoidIntersectionImageFilter()
        rei.SetDensity(self.get_density())
        rei.SetAngle(self.get_angle())
        rei.SetCenter(self.get_center(t))  #
        rei.SetAxis(self.get_axis())

        rei.SetGeometry(geometry)
        reiImage = rei.Execute(empty_image_detector)

        return srtk.GetArrayFromImage(reiImage)[0, 0, :]

class AnalyticalEllipse(object):
    """
        Class of ellipse where projections are computed according to analytical (i.e. exact)
        formulas. Hence, this is not a simulation but a computation.
    """

    def __init__(self, params):
        self.params = params

    def get_density(self):
        return self.params.ellipseDensity

    def get_angle(self):
        return self.params.ellipseAngle

    def get_center(self, t):
        """
            Since it is moving, the position depends on t
        """
        T = self.params.T
        v = self.params.v
        v2 = self.params.v2

        return [self.params.ellipseCenterX + (t + T / 2) * v,
                self.params.ellipseCenterY + (t + T / 2) * v2,
                0]

    def get_axis(self):
        return [self.params.ellipseSemiAxisX,
                self.params.ellipseSemiAxisY,
                self.params.ellipseSemiAxisY]

    def compute_projection(self, t, source, detector):
        """
            Simulate fan-beam acquisition of the object with given
            source and detector, at time t
        """
        N = self.params.imageSize
        results = np.zeros(N)

        # general parameters
        alpha = self.params.get_alpha_range()
        omega = self.params.omega / 360 * 2 * np.pi

        # ellipse parameters
        x, y, _ = self.get_center(t)
        a, b, _ = self.get_axis()
        if (a != b):
            raise ValueError("Ellipse is not a circle (the analytical formula only works with circle)", a, b)

        s1, s2 = source.get_position(t)

        # for i in np.arange(N):
        # 	# i is the number of the pixel in the image printed on the detector
        # 	# there is no "resolution" parameter, meaning that there is 1 pixel
        # 	# per millimeter
        # 	# TODO : add this parameter?

        # 	# phi is the angle between the beam and the y-axis
        # 	phi = omega*t + alpha[i]

        # 	# computation of the distance between the center of the circle
        # 	# and the beam
        # 	dist = np.abs( (s1-x)*np.cos(phi) + (s2-y)*np.sin(phi) )

        # 	# stores in the array
        # 	if dist > a:
        # 		results[i] = 0
        # 	else:
        # 		results[i] = 2 * np.sqrt(a**2 - dist**2)

        # phi is the angle between the beam and the y-axis
        phi = omega * t + alpha

        # distance between the center of the circle and the beam
        dist = np.abs((s1 - x) * np.cos(phi) + (s2 - y) * np.sin(phi))

        # results = (dist<a) * 2 * np.sqrt(a**2 - dist**2)
        # [x if x < 5 else 0 for x in np.arange(10)]
        # ipdb.set_trace()
        results[dist < a] = (2 * np.sqrt(a ** 2 - dist ** 2))[dist < a]

        return self.get_density() * results