import SimpleRTK as srtk

class MovingEllipse(object):
    """docstring for MovingEllipse"""

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
