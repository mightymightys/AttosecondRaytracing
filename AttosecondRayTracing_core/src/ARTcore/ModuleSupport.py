"""
Provides classes for different shapes of a support for optics.

A support fixes the spatial extent and  outer shape of optics - think of it as a substrate.
Technically, the support is a section of the x-y-plane in the element coordinate frame.

The support will become an attribute of Mirror or Mask objects.

![Illustration of different kinds of support.](supports.svg)



Created in Sept 2020

@author: Anthony Guillaume and Stefan Haessler
"""

# %%

import numpy as np
import math
import ARTcore.ModuleGeometry as mgeo
from abc import ABC, abstractmethod
import logging

logger = logging.getLogger(__name__)


# %% Abstract Support Class
class Support(ABC):
    """
    Abstract base class for optics supports.
    Supports the `in` operator to test points.
    """
    @abstractmethod
    def _sdf(self, Point):
        """
        Signed distance function for the support.
        """
        pass

    def __contains__(self, Point):
        """
        Interface allowing the use of the `in` operator to check if a point is within the support.
        """
        return self._sdf(Point) <= 0
    
    def _estimate_size(self, initial_distance=1000):
        """
        Estimate the size of the support using the signed distance function.
        """
        directions = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])  # Simple axis-aligned directions
        points = directions * initial_distance  # Generate points on a circle of the current radius
        distances = np.array([self._sdf(p) for p in points])
        avg_distance = np.mean(distances)

        return initial_distance - avg_distance
    


    

# %% Round Support
class SupportRound(Support):
    """
    A round support for optics.

    Attributes
    ----------
        radius : float
            The radius of the support in mm.
    """

    def __init__(self, Radius: float):
        """
        Create a round support.

        Parameters
        ----------
            Radius : float
                The radius of the support in mm.

        """
        self.radius = Radius

    def __repr__(self):
        return f"{type(self).__name__}({self.radius})"

    def _sdf(self, Point):
        """
        Return signed distance from the support.
        """
        return mgeo.SDF_Circle(Point, self.radius)
    
    def _get_grid(self, NbPoint: int, **kwargs):
        """
        Return a Point array with the coordinates of a number NbPoints of points.
        The points are distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        return mgeo.Point(np.hstack((MatrixXY, np.zeros((MatrixXY.shape[0], 1)))))
    
    def _get_edges(self, NbPoint=100):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        return flatten_point_arrays(gen_circle_contour(self.radius, NbPoint), [])


# %% Round Support with Hole
class SupportRoundHole(Support):
    """
    A round support for optics with a round hole.

    Attributes
    ----------
        radius : float
            The radius of the support in mm.

        radiushole : float
            The radius of the hole in mm.

        centerholeX : float
            The x-cordinate of the hole's centre, in mm.

        centerholeY : float
            The y-cordinate of the hole's centre, in mm.

    """

    def __init__(self, Radius: float, RadiusHole: float, CenterHoleX: float, CenterHoleY: float):
        """
        Parameters
        ----------
        Radius : float
            The radius of the support in mm.

        RadiusHole : float
            The radius of the hole in mm.

        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.

        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.

        """
        self.radius = Radius
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY

    def __repr__(self):
        return f"{type(self).__name__}(Radius={self.radius}, RadiusHole={self.radiushole}, CenterHoleX={self.centerholeX}, CenterHoleY={self.centerholeY})"

    def _sdf(self, Point):
        """Return signed distance from the support."""
        support = mgeo.SDF_Circle(Point, self.radius)
        hole =  mgeo.SDF_Circle(Point[:2] - np.array([self.centerholeX, self.centerholeY]), self.radiushole)
        return mgeo.Difference_SDF(support, hole)
    
    def _get_grid(self, NbPoint, **kwargs):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        return MatrixXY
        

    def _ContourSupport(self, Figure):
        """Draws support contour in MirrorProjection plots."""
        raise ValueError("No visualisation")

    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        N_outer = int(round(NbPoint - NbPoint * self.radiushole / self.radius))
        outer = gen_circle_contour(self.radius, N_outer)
        hole = gen_circle_contour(self.radiushole, NbPoint - N_outer) + np.array([self.centerholeX, self.centerholeY])
        return flatten_point_arrays(outer, [hole], edges=edges)


# %% Rectangular Support
class SupportRectangle(Support):
    """
    A rectangular support for optics.

    Attributes
    ----------
        dimX : float
            The dimension in mm along x.

        dimY : float
            The dimension in mm along y.

    """

    def __init__(self, DimensionX: float, DimensionY: float):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.

        DimensionY : float
            The dimension in mm along y.

        """
        self.dimX = DimensionX
        self.dimY = DimensionY

    def __repr__(self):
        return f"{type(self).__name__}(DimensionX={self.dimX}, DimensionY={self.dimY})"

    def _sdf(self, Point):
        """Return signed distance from the support."""
        return mgeo.SDF_Rectangle(Point, self.dimX, self.dimY)

    def _get_grid(self, NbPoints: int, **kwargs) -> list[np.ndarray]:
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        nbx = int(
            np.sqrt(self.dimX / self.dimY * NbPoints + 0.25 * (self.dimX - self.dimY) ** 2 / self.dimY**2)
            - 0.5 * (self.dimX - self.dimY) / self.dimY
        )
        nby = int(NbPoints / nbx)
        x = np.linspace(-self.dimX / 2, self.dimX / 2, nbx)
        y = np.linspace(-self.dimY / 2, self.dimY / 2, nby)
        ListCoordXY = []
        for i in x:
            for j in y:
                ListCoordXY.append(np.array([i, j]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        """Draws support contour in MirrorProjection plots."""
        raise ValueError("No visualisation")


    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        return flatten_point_arrays(gen_rectangle_contour(self.dimX, self.dimY, NbPoint), [], edges=edges)


# %% Rectangular Support with Hole
class SupportRectangleHole(Support):
    """
    A rectangular support for optics with a round hole.

    Attributes
    ----------
        dimX : float
            The dimension in mm along x.

        dimY : float
            The dimension in mm along y.

        radiushole : float
            The radius of the hole in mm.

        centerholeX : float
            The x-cordinate of the hole's centre, in mm.

        centerholeY : float
            The y-cordinate of the hole's centre, in mm.

    """

    def __init__(self, DimensionX: float, DimensionY: float, RadiusHole: float, CenterHoleX: float, CenterHoleY: float):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.

        DimensionY : float
            The dimension in mm along y.

        RadiusHole : float
            The radius of the hole in mm.

        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.

        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.

        """
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.radiushole = RadiusHole
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY

    def __repr__(self):
        return f"{type(self).__name__}(DimensionX={self.dimX}, DimensionY={self.dimY}, RadiusHole={self.radiushole}, CenterHoleX={self.centerholeX}, CenterHoleY={self.centerholeY})"

    def _sdf(self, Point):
        """Return signed distance from the support."""
        support = mgeo.SDF_Rectangle(Point, self.dimX, self.dimY)
        hole = mgeo.SDF_Circle(Point[:2] - np.array([self.centerholeX, self.centerholeY]), self.radiushole)
        return mgeo.Difference_SDF(support, hole)
    


# %% Rectangular Support with Rectangular Hole
class SupportRectangleRectHole(Support):
    """
    A rectangular support for optics, with a rectangular hole.

    Attributes
    ----------
        dimX : float
            The dimension in mm along x.

        dimY : float
            The dimension in mm along y.

        holeX : float
            The dimension of the hole in mm along x.

        holeY : float
            The dimension of the hole in mm along y.

        centerholeX : float
            The x-cordinate of the hole's centre, in mm.

        centerholeY : float
            The y-cordinate of the hole's centre, in mm.

    """

    def __init__(
        self, DimensionX: float, DimensionY: float, HoleX: float, HoleY: float, CenterHoleX: float, CenterHoleY: float
    ):
        """
        Parameters
        ----------
        DimensionX : float
            The dimension in mm along x.

        DimensionY : float
            The dimension in mm along y.

        HoleX : float
            The dimension of the hole in mm along x.

        HoleY : float
            The dimension of the hole in mm along y.

        CenterHoleX : float
            The x-cordinate of the hole's centre, in mm.

        CenterHoleY : float
            The y-cordinate of the hole's centre, in mm.

        """
        self.dimX = DimensionX
        self.dimY = DimensionY
        self.holeX = HoleX
        self.holeY = HoleY
        self.centerholeX = CenterHoleX
        self.centerholeY = CenterHoleY

    def __repr__(self):
        return f"{type(self).__name__}(DimensionX={self.dimX}, DimensionY={self.dimY}, HoleX={self.holeX}, HoleY={self.holeY}, CenterHoleX={self.centerholeX}, CenterHoleY={self.centerholeY})"

    def _sdf(self, Point):
        """Return signed distance from the support."""
        support = mgeo.SDF_Rectangle(Point, self.dimX, self.dimY)
        hole = mgeo.SDF_Rectangle(Point[:2] - np.array([self.centerholeX, self.centerholeY]), self.holeX, self.holeY)
        return mgeo.Difference_SDF(support, hole)
