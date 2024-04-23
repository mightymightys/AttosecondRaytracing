"""
Provides a class for masks, which are a type of optical element that simply stops rays that hit it.

Also provides the function *TransmitMaskRayList* that returns the rays transmitted by the mask.

![Illustration the Mask-class.](Mask.svg)



Created in 2021

@author: Stefan Haessler
"""
# %% Modules

import numpy as np
import ARTcore.ModuleGeometry as mgeo


# %%############################################################################
class Mask:
    """
    A mask: a plane surface of the shape of its [Support](ModuleSupport.html), which stops all rays that hit it,
    while all other rays continue their path unchanged.

    Attributes
    ----------
        support : [Support](ModuleSupport.html)-object

        type : str 'Ellipsoidal Mirror'.

    Methods
    ----------
        get_normal(Point)

        get_centre()

        get_grid3D(NbPoints)

    """

    def __init__(self, Support):
        """
        Parameters
        ----------
            Support : [Support](ModuleSupport.html)-object
        """
        self.type = "Mask"
        self.support = Support

    def _get_intersection(self, Ray):
        """Return the intersection point between Ray and the xy-plane,
        where rays do NOT hit the mask's support"""
        t = -Ray.point[2] / Ray.vector[2]
        I = Ray.vector * t + Ray.point
        if t > 0 and not self.support._IncludeSupport(I):
            PointIntersection = I
        else:
            PointIntersection = None

        return PointIntersection

    def get_normal(self, Point):
        """Return normal unit vector in point 'Point' on the mask."""
        Normal = np.array([0, 0, 1])
        return Normal

    def get_centre(self):
        """Return 3D coordinates of the point on the mask surface at the center of its support."""
        return np.array([0, 0, 0])

    def get_grid3D(self, NbPoint, **kwargs):
        """
        Returns list of numpy-arrays containing the 3D-coordinates of points in the mirror surface,
        sampling the support in a number NbPoints of points.
        """
        E = "edges" in kwargs and kwargs["edges"]
        ListCoordXYZ = []
        contour = int(round(0.1 * NbPoint))
        contours = self.support._Contour_points(contour, **kwargs)
        if E:
            contours, contour_edges = contours
        ListCoordXY = contours + self.support._get_grid(NbPoint - contour)
        for k in ListCoordXY:
            z = 0
            ListCoordXYZ.append(np.array([k[0], k[1], z]))
        if E:
            return ListCoordXYZ, contour_edges
        return ListCoordXYZ


# %%############################################################################
def _TransmitMaskRay(Mask, PointMask, Ray):
    """Returns the transmitted ray"""
    PointRay = Ray.point
    VectorRay = Ray.vector
    NormalMask = Mask.get_normal(PointMask)

    AngleIncidence = mgeo.AngleBetweenTwoVectors(VectorRay, NormalMask)
    Path = Ray.path + (np.linalg.norm(PointMask - PointRay),)

    RayTransmitted = Ray.copy_ray()
    RayTransmitted.point = PointMask
    RayTransmitted.vector = VectorRay
    RayTransmitted.incidence = AngleIncidence
    RayTransmitted.path = Path

    return RayTransmitted


# %%
def TransmitMaskRayList(Mask, RayList):
    """
    Returns the the transmitted rays that pass the mask for the list of
    incident rays ListRay.

    Rays that hit the support are not further propagated.

    Updates the reflected rays' incidence angle and path.

    Parameters
    ----------
        Mask : Mask-object

        ListRay : list[Ray-object]

    """
    ListRayTransmitted = []
    for Ray in RayList:
        PointMask = Mask._get_intersection(Ray)

        if PointMask is not None:
            RayTransmitted = _TransmitMaskRay(Mask, PointMask, Ray)
            ListRayTransmitted.append(RayTransmitted)

    return ListRayTransmitted
