"""
Provides classes for different reflective surfaces.

Created in Nov 2024

@author: Andre Kalouguine
"""
import ARTcore.ModuleGeometry as mgeo

import numpy as np
from copy import copy
from abc import ABC, abstractmethod


class Surface(ABC):
    """
    Abstract base class for surfaces.
    
    This is where we should ideally define the roughness, scratches/digs, spectral response, etc.

    As for the global sag, that might be better defined in the optical elements themselves.
    """
    pass

class IdealSurface(Surface):
    """
    Ideal surface, i.e. no roughness, no scratches, no digs, no spectral response.
    """
    def __init__(self):
        super().__init__()
    def reflect_ray(self, ray, point, normal):
        VectorRay = ray.vector
        VectorRayReflected = VectorRay-2*normal*np.dot(VectorRay,normal) # Is it any better than SymmetricalVector?
        local_reflectedray = copy(ray)
        local_reflectedray.point = point
        local_reflectedray.vector = VectorRayReflected
        local_reflectedray.incidence = mgeo.AngleBetweenTwoVectors(-VectorRay, normal)
        local_reflectedray.path = ray.path + (np.linalg.norm(point - ray.point),)
        return local_reflectedray
    