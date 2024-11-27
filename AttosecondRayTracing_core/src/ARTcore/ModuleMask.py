"""
Provides a class for masks, which are a type of optical element that simply stops rays that hit it.

Also provides the function *TransmitMaskRayList* that returns the rays transmitted by the mask.

![Illustration the Mask-class.](Mask.svg)



Created in 2021

@author: Stefan Haessler + André Kalouguine
"""
# %% Modules
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin
import ARTcore.ModuleSupport as msup
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleMirror as mmirror

import numpy as np
from copy import copy
import logging

logger = logging.getLogger(__name__)

# %%############################################################################
class Mask(moe.OpticalElement):
    """
    A mask: a plane surface of the shape of its [Support](ModuleSupport.html), which stops all rays that hit it,
    while all other rays continue their path unchanged. No diffraction effects.

    Attributes
    ----------
        support : [Support](ModuleSupport.html)-object

        type : str 'Mask'.

    Methods
    ----------
        get_local_normal(Point)

        get_centre()

        get_grid3D(NbPoints)

    """

    def __init__(self, Support, **kwargs):
        """
        Parameters
        ----------
            Support : [Support](ModuleSupport.html)-object
        """
        self.type = "Mask"
        self.support = Support
        self.curvature = mmirror.Curvature.FLAT

        self.r0 = mgeo.Point([0.0, 0.0, 0.0])
        self._r = mgeo.Vector([0.0, 0.0, 0.0])
        self._q = np.quaternion(1)

        self.centre_ref = mgeo.Point([0, 0, 0])

        self.support_normal_ref = mgeo.Vector([0, 0, 1])
        self.normal_ref = mgeo.Vector([0, 0, 1])
        self.majoraxis_ref = mgeo.Vector([1, 0, 0])

        self.add_global_vectors("support_normal", "normal", "majoraxis")
        self.add_global_points("centre")

        super().__init__()

    def _get_intersection(self, Ray):
        """
        Return the intersection point between Ray and the xy-plane.
        If not in alignment mode, any intersection point that is on the support is blocked.
        """
        t = -Ray.point[2] / Ray.vector[2]
        I = mgeo.Point(Ray.point +Ray.vector * t)
        return I, I in self.support
    
    def _get_intersections(self, RayList):
        """
        Return the intersection point between Ray and the xy-plane.
        If not in alignment mode, any intersection point that is on the support is blocked.
        """
        vector = mgeo.VectorArray([ray.vector for ray in RayList])
        point = mgeo.PointArray([ray.point for ray in RayList])
        t = -point[:, 2] / vector[:, 2]
        I = mgeo.PointArray(point + vector * t[:, np.newaxis])
        return I, [not(i in self.support) for i in I]

    def get_local_normal(self, Point):
        """Return normal unit vector in point 'Point' on the mask."""
        return np.array([0, 0, 1])
    
    def _zfunc(self, PointArray):
        """Return the z-values of the plane surface at the points in PointArray."""
        return np.zeros(len(PointArray))

    def propagate_raylist(self, RayList, alignment=False):
        """
        Propagate a list of rays through the mask, returning the transmitted rays.

        Parameters
        ----------
        RayList : list
            List of Ray objects to be propagated through the mask.
        alignment : bool, optional
            If True, the alignment of the mask is considered. Default is False.

        Returns
        -------
        list
            List of transmitted Ray objects.
        """
        transmitted_rays = []
        local_rays = [ray.to_basis(*self.basis) for ray in RayList]
        intersection_points, OK = self._get_intersections(local_rays)
        N = len(RayList)
        for i in range(N):
            # Then we find the intersection point in the mirror reference frame
            intersection_point = intersection_points[i]
            local_ray = local_rays[i]
            if OK[i]:
                local_transmittedray = copy(local_ray)
                local_transmittedray.point = intersection_point
                local_transmittedray.vector = local_ray.vector
                local_transmittedray.incidence = mgeo.AngleBetweenTwoVectors(-local_ray.vector, self.normal_ref)
                local_transmittedray.path = local_ray.path + (np.linalg.norm(intersection_point - local_ray.point),)
                transmitted_rays.append(local_transmittedray.from_basis(*self.basis))
        if len(transmitted_rays) == 0:
            logger.warning("No rays were transmitted by the mask.")
            logger.debug(f"Mask: {self}")
            logger.debug(f"First ray: {RayList[0]}")
            logger.debug(f"First ray in mask reference frame: {RayList[0].to_basis(self.r0, self.r, self.q)}")
            logger.debug(f"First ray intersection point: {self._get_intersection(RayList[0].to_basis(self.r0, self.r, self.q))}")
        return mray.RayList.from_list(transmitted_rays)


# %%############################################################################
def _TransmitMaskRay(Mask, PointMask, Ray):
    """Returns the transmitted ray"""
    PointRay = Ray.point
    VectorRay = Ray.vector
    NormalMask = Mask.get_local_normal(PointMask)

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
