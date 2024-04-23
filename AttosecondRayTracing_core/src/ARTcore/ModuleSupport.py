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


# %%


class Support(ABC):
    """Abstract base class for optics supports."""

    @abstractmethod
    def _IncludeSupport(self, Point):
        pass

    @abstractmethod
    def _get_grid(self, NbPoint, **kwargs):
        pass

    @abstractmethod
    def _ContourSupport(self, Figure):
        pass


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

    def _IncludeSupport(self, Point):
        """Return whether 2D-Point is within the rectangular Support."""
        return mgeo.IncludeDisk(self.radius, Point)

    def _get_grid(self, NbPoint: int, **kwargs):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k, 0]
            y = MatrixXY[k, 1]
            ListCoordXY.append(np.array([x, y]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        raise ValueError("No visualisation")

    def _CircumRect(self):
        return np.array([self.radius * 2, self.radius * 2])

    def _CircumCirc(self):
        return self.radius

    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        return flatten_point_arrays(gen_circle_contour(self.radius, NbPoint), [], edges=edges)


# %%
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

    def _IncludeSupport(self, Point):
        """Returns whether 2D-Point is within the rectangular Support."""
        return mgeo.IncludeDisk(self.radius, Point) and not (
            mgeo.IncludeDisk(self.radiushole, Point - np.array([self.centerholeX, self.centerholeY, 0]))
        )

    def _get_grid(self, NbPoint, **kwargs):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a Vogel-spiral on the support, with the origin in the center of the support.
        """
        MatrixXY = mgeo.SpiralVogel(NbPoint, self.radius)
        ListCoordXY = []
        for k in range(NbPoint):
            x = MatrixXY[k, 0]
            y = MatrixXY[k, 1]
            if (x - self.centerholeX) ** 2 + (y - self.centerholeY) ** 2 > self.radiushole**2:
                ListCoordXY.append(np.array([x, y]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        """Draws support contour in MirrorProjection plots."""
        raise ValueError("No visualisation")

    def _CircumRect(self):
        return np.array([self.radius * 2, self.radius * 2])

    def _CircumCirc(self):
        return self.radius

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


# %%


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

    def _IncludeSupport(self, Point: np.ndarray) -> bool:
        """Return whether 2D-Point is within the rectangular Support."""
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point)

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

    def _CircumRect(self):
        return np.array([self.dimX, self.dimY])

    def _CircumCirc(self):
        return np.sqrt(self.dimX**2 + self.dimY**2) / 2

    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        return flatten_point_arrays(gen_rectangle_contour(self.dimX, self.dimY, NbPoint), [], edges=edges)


# %%
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

    def _IncludeSupport(self, Point):
        """Returns whether 2D-Point is within the rectangular Support."""
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not (
            mgeo.IncludeDisk(self.radiushole, Point - np.array([self.centerholeX, self.centerholeY, 0]))
        )

    def _get_grid(self, NbPoint, **kwargs):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        x = np.linspace(-self.dimX / 2, self.dimX / 2, int(self.dimX / self.dimY * np.sqrt(NbPoint)))
        y = np.linspace(-self.dimY / 2, self.dimY / 2, int(self.dimY / self.dimX * np.sqrt(NbPoint)))

        ListCoordXY = []
        for i in x:
            for j in y:
                if (i - self.centerholeX) ** 2 + (j - self.centerholeY) ** 2 > self.radiushole**2:
                    ListCoordXY.append(np.array([i, j]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        """Draws support contour in MirrorProjection plots."""
        raise ValueError("No visualisation")

    def _CircumRect(self):
        return np.array([self.dimX, self.dimY])

    def _CircumCirc(self):
        return np.sqrt(self.dimX**2 + self.dimY**2) / 2

    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        outer_length = 2 * (self.dimX + self.dimY)
        hole_length = 2 * np.pi * self.radiushole
        total_length = outer_length + hole_length
        NbHole = int(round(hole_length / total_length * NbPoint))
        outer = gen_rectangle_contour(self.dimX, self.dimY, NbPoint - NbHole)
        hole = gen_circle_contour(self.radiushole, NbHole) + np.array([self.centerholeX, self.centerholeY])
        return flatten_point_arrays(outer, [hole], edges=edges)


# %%
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

    def _IncludeSupport(self, Point):
        """Return whether 2D-Point is within the rectangular Support."""
        return mgeo.IncludeRectangle(self.dimX, self.dimY, Point) and not mgeo.IncludeRectangle(
            self.holeX, self.holeY, Point - np.array([self.centerholeX, self.centerholeY, 0])
        )

    def _get_grid(self, NbPoint, **kwargs):
        """
        Returns a list of 2D-numpy-arrays with the coordinates a number NbPoints of points,
        distributed as a regular grid on the support, with the origin in the center of the support.
        """
        nbx = int(
            np.sqrt(self.dimX / self.dimY * NbPoint + 0.25 * (self.dimX - self.dimY) ** 2 / self.dimY**2)
            - 0.5 * (self.dimX - self.dimY) / self.dimY
        )
        nby = int(NbPoint / nbx)
        x = np.linspace(-self.dimX / 2, self.dimX / 2, nbx)
        y = np.linspace(-self.dimY / 2, self.dimY / 2, nby)

        ListCoordXY = []
        for i in x:
            for j in y:
                if abs(i - self.centerholeX) > self.holeX / 2 or abs(j - self.centerholeY) > self.holeY / 2:
                    ListCoordXY.append(np.array([i, j]))
        return ListCoordXY

    def _ContourSupport(self, Figure):
        """Draws support contour in MirrorProjection plots."""
        raise ValueError("No visualisation")

    def _CircumRect(self):
        return np.array([self.dimX, self.dimY])

    def _CircumCirc(self):
        return np.sqrt(self.dimX**2 + self.dimY**2) / 2

    def _Contour_points(self, NbPoint=100, edges=False):
        """
        Return a list of 2D-numpy-arrays with the coordinates a number NbPoints of points.

        The points are distributed along the contour of the support so as to clearly define the edges.
        The putpose is to use these points to draw a nicer mesh of the mirrors.
        """
        outer_length = 2 * (self.dimX + self.dimY)
        hole_length = 2 * (self.holeX + self.holeY)
        total_length = outer_length + hole_length
        NbHole = int(round(hole_length / total_length * NbPoint))
        outer = gen_rectangle_contour(self.dimX, self.dimY, NbPoint - NbHole)
        hole = gen_rectangle_contour(self.holeX, self.holeY, NbHole) + np.array([self.centerholeX, self.centerholeY])
        return flatten_point_arrays(outer, [hole[::-1]], edges=edges)


def find_hull(points):
    # start from leftmost point
    current_point = min(range(len(points)), key=lambda i: points[i][0])
    # initialize hull with current point
    hull = [current_point]
    # initialize list of linked points
    linked = []
    # continue until all points have been linked
    while len(linked) < len(points) - 1:
        # initialize minimum distance and closest point
        min_distance = math.inf
        closest_point = None
        # find closest unlinked point to current point
        for i, point in enumerate(points):
            if i not in linked:
                distance = math.dist(points[current_point], point)
                if distance < min_distance:
                    min_distance = distance
                    closest_point = i
        # add closest point to hull and linked list
        hull.append(closest_point)
        linked.append(closest_point)
        # update current point
        current_point = closest_point
    # add link between last point and first point
    hull.append(hull[0])
    # convert hull to a list of pairs of indices
    indices = [[hull[i], hull[i + 1]] for i in range(len(hull) - 1)]
    return indices


def gen_rectangle_contour(dimX, dimY, NbPoints):
    # calculate number of points on each side
    nX = math.ceil(dimX / (dimX + dimY) * NbPoints)
    nY = NbPoints - nX
    # calculate distance between points on each side
    dX = dimX / (nX - 1)
    dY = dimY / (nY - 1)
    # generate points on top side
    top_points = [(i * dX - dimX / 2, dimY / 2) for i in range(nX)]
    # generate points on right side
    right_points = [(dimX / 2, dimY / 2 - i * dY) for i in range(1, nY)]
    # generate points on bottom side
    bottom_points = [(dimX / 2 - i * dX, -dimY / 2) for i in range(1, nX)]
    # generate points on left side
    left_points = [(-dimX / 2, -dimY / 2 + i * dY) for i in range(1, nY - 1)]
    # concatenate points in counterclockwise order
    points = top_points + right_points + bottom_points + left_points
    return np.array(points)


def gen_circle_contour(radius, NbPoints):
    # calculate angle between points
    if NbPoints == 0:
        return np.array([])
    else:    
        angle = 2 * math.pi / NbPoints
    # generate points
    points = [(radius * math.cos(i * angle), radius * math.sin(i * angle)) for i in range(NbPoints)]
    return np.array(points)


def flatten_point_arrays(outer, holes=[], edges=False):
    coords = [i for i in outer]
    edges_list = []
    offset = len(coords)
    for arr in holes:
        # add coordinates to list
        coords.extend(arr)
        # add edges to list
        if edges:
            n = arr.shape[0]
            edge_indices = np.arange(n)
            edge_indices += offset
            edges_list += [edge_indices.T.tolist() + [offset]]
        # update offset
        offset += arr.shape[0]
    if edges:
        return coords, edges_list
    else:
        return coords
