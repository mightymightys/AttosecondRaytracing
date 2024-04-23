"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
# %% Modules

import numpy as np

# %%
class Ray:
    """
    Represents an optical ray from geometrical optics.
    """

    # __slots__ = ('_point', '_vector', '_path', '_number', '_wavelength', '_incidence', '_intensity') #doesn't bring a significant performance improvement

    def __init__(
        self,
        Point: np.ndarray,
        Vector: np.ndarray,
        Path=(0.0,),
        Number=None,
        Wavelength=None,
        Incidence=None,
        Intensity=None,
    ):
        """
        Parameters
        ----------
            Point : np.ndarray
                Point source of the ray..

            Vector : np.ndarray
                Direction unit-vector of the ray..

            Path : List of floats, optional TODO correct considering update to list
                Path length in mm covered by ray *before* its starting Point. The default is 0.

            Number : int, optional
                Index needed to link the incident ray and the reflected ray. The default is None, but it should be set by afterwards setting the attribute Ray._number.!

            Wavelength : float, optional
                Wavelength of the light in mm. The default is None.

            Incidence : float, optional
                Incidence angle in radian of the ray on the optical element *from which it originates*.
                Gets calculated during ray tracing. The default is None. The source rays therefore have "None" incidence.

            Intensity : float, optional
                Intensity, or more strictly speaking the fluence fraction in arb.u. carried by the ray. The default is None.

        """
        self.point = Point  # the setter checks
        self.vector = Vector  # the setter checks and normalizes
        self._path = Path
        self._wavelength = Wavelength
        self._incidence = Incidence
        self._intensity = Intensity
        if type(Number) == int or Number is None:
            self._number = Number
        else:
            raise TypeError("Ray Number must be an integer.")

    """
    Removing the setters and getters with their checks would gain about 10-30% speed because we use the Ray-class a lot.
    Then the main risk seems that we can no longer trust in the vector being a unit vector.
    """

    @property
    def point(self):
        return self._point

    @point.setter
    def point(self, Point):
        if type(Point) == np.ndarray and len(Point) == 3:
            self._point = Point
        else:
            raise TypeError("Ray Point must be a 3D numpy.ndarray, but it is  %s." % type(Point))

    @property
    def vector(self):
        return self._vector

    @vector.setter
    def vector(self, Vector):
        if type(Vector) == np.ndarray and len(Vector) == 3 and np.linalg.norm(Vector) > 1e-9:
            self._vector = Vector / np.linalg.norm(Vector)
        else:
            raise TypeError("Ray Vector must be a 3D numpy.ndarray with finite length.")

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, Path):
        self._path = Path

    @property
    def number(self):
        return self._number

    # actually, we don't want to number to be changable afterwards, so not setter
    # @number.setter
    # def number(self, Number):
    #     if type(Number) == int or Number is None:
    #           self._number = Number
    #     else: raise TypeError('Ray Number must be an integer.')

    @property
    def wavelength(self):
        return self._wavelength

    @wavelength.setter
    def wavelength(self, Wavelength):
        if type(Wavelength) in [int, float, np.float64]:
            self._wavelength = Wavelength
        else:
            raise TypeError("Ray Wavelength must be int or float or None.")

    @property
    def incidence(self):
        return self._incidence

    @incidence.setter
    def incidence(self, Incidence):
        if type(Incidence) in [float, np.float64]:
            self._incidence = Incidence
        else:
            raise TypeError("Ray Incidence must be a float or None.")

    @property
    def intensity(self):
        return self._intensity

    @intensity.setter
    def intensity(self, Intensity):
        if type(Intensity) in [int, float, np.float64]:
            self._intensity = Intensity
        else:
            raise TypeError("Ray Intensity must be int or float or None.")

    # %%
    def copy_ray(self):
        """
        Returns a new OpticalRay object with the same properties.
        """
        return Ray(self.point, self.vector, self.path, self.number, self.wavelength, self.incidence, self.intensity)

    def __hash__(self):
        point_tuple = tuple(self.point.reshape(1, -1)[0])
        vector_tuple = tuple(self.vector.reshape(1, -1)[0])
        return hash(
            point_tuple + vector_tuple + (self.path, self.number, self.wavelength, self.incidence, self.intensity)
        )
