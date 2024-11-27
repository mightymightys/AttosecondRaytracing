"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler + André Kalouguine
"""
# %% Modules

import numpy as np
import logging
import ARTcore.ModuleGeometry as mgeo
logger = logging.getLogger(__name__)
import time
from dataclasses import dataclass, replace

# %% Ray class definition

@dataclass(slots=True)
class Ray:
    """
    Represents an optical ray from geometrical optics.
    In order to optimize lookups, Rays are represented as dataclasses.
    The mandatory attributes are:
    - point: mgeo.Point
    - vector: mgeo.Vector
    """
    point: mgeo.Point
    vector: mgeo.Vector
    path: tuple = (0.0,)
    number: int = None
    wavelength: float = None
    incidence: float = None
    intensity: float = None

    def to_basis(self, r0, r, q):
        """
        Transforms the ray to the basis defined by the origin r0, the x-axis r and the y-axis q.
        """
        return Ray(
            self.point.to_basis(r0, r, q),
            self.vector.to_basis(r0, r, q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )
    
    def from_basis(self, r0, r, q):
        """
        Transforms the ray from the basis defined by the origin r0, the x-axis r and the y-axis q.
        """
        return Ray(
            self.point.from_basis(r0, r, q),
            self.vector.from_basis(r0, r, q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )

    def rotate(self, q):
        """
        Rotates the ray by the quaternion q.
        """
        return Ray(
            self.point.rotate(q),
            self.vector.rotate(q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )
    
    def translate(self, t):
        """
        Rotates the ray by the vector t
        """
        return Ray(
            self.point.translate(t),
            self.vector,
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )

    def __copy__(self):
        """
        Returns a new OpticalRay object with the same properties.
        """
        return replace(self)

    def __hash__(self):
        point_tuple = tuple(self.point.reshape(1, -1)[0])
        vector_tuple = tuple(self.vector.reshape(1, -1)[0])
        return hash(
            point_tuple + vector_tuple + (self.path, self.number, self.wavelength, self.incidence, self.intensity)
        )


@dataclass(slots=True)
class RayList:
    """
    Class representing a list of rays in a form that is well suited to calculations.
    Specifically, the points, vectors etc are all numpy arrays
    It can then be converted to a list of Ray objects with the method 
    """
    point: mgeo.PointArray
    vector: mgeo.VectorArray
    path: np.ndarray
    number: np.ndarray
    wavelength: np.ndarray
    incidence: np.ndarray
    intensity: np.ndarray

    @classmethod
    def from_list(cls, rays):
        N = len(rays)
        N_r = len(rays[0].path)
        points = np.zeros((N, 3))
        vectors = np.zeros((N, 3))
        paths = np.zeros((N, N_r))
        numbers = np.zeros(N, dtype=int)
        wavelengths = np.zeros(N)
        incidences = np.zeros(N)
        intensities = np.zeros(N)
        for i, ray in enumerate(rays):
            points[i:,] = ray.point
            vectors[i:,] = ray.vector
            paths[i:,] = ray.path
            numbers[i] = ray.number
            wavelengths[i] = ray.wavelength
            incidences[i] = ray.incidence
            intensities[i] = ray.intensity
        return cls(
            mgeo.PointArray(points),
            mgeo.VectorArray(vectors),
            paths,
            numbers,
            wavelengths,
            incidences,
            intensities,
        )

    @property
    def N(self):
        return len(self.point)
    
    def __getitem__(self, key):
        if isinstance(key, slice) or isinstance(key, list):
            return RayList(
                self.point[key],
                self.vector[key],
                self.path[key],
                self.number[key],
                self.wavelength[key],
                self.incidence[key],
                self.intensity[key],
            )
        return Ray(
            self.point[key],
            self.vector[key],
            self.path[key],
            self.number[key],
            self.wavelength[key],
            self.incidence[key],
            self.intensity[key],
        )
    
    def __len__(self):
        return self.N
    
    def __iter__(self):
        return (Ray(point=self.point[i],
                    vector=self.vector[i],
                    path=tuple(self.path[i]),
                    number=int(self.number[i]),
                    wavelength=float(self.wavelength[i]),
                    incidence=float(self.incidence[i]),
                    intensity=float(self.intensity[i])) for i in range(self.N))

    def to_basis(self, r0, r, q):
        """
        Transforms the ray to the basis defined by the origin r0, the x-axis r and the y-axis q.
        """
        return RayList(
            self.point.to_basis(r0, r, q),
            self.vector.to_basis(r0, r, q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )
    
    def from_basis(self, r0, r, q):
        """
        Transforms the ray from the basis defined by the origin r0, the x-axis r and the y-axis q.
        """
        return RayList(
            self.point.from_basis(r0, r, q),
            self.vector.from_basis(r0, r, q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )

    def rotate(self, q):
        """
        Rotates the ray by the quaternion q.
        """
        return Ray(
            self.point.rotate(q),
            self.vector.rotate(q),
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )
    
    def translate(self, t):
        """
        Rotates the ray by the vector t
        """
        return Ray(
            self.point.translate(t),
            self.vector,
            self.path,
            self.number,
            self.wavelength,
            self.incidence,
            self.intensity,
        )

    def __copy__(self):
        """
        Returns a new OpticalRay object with the same properties.
        """
        return RayList(
            self.point.copy(),
            self.vector.copy(),
            self.path.copy(),
            self.number.copy(),
            self.wavelength.copy(),
            self.incidence.copy(),
            self.intensity.copy(),
        )

    def __hash__(self):
        points = self.point.reshape(1, -1)[0]
        vectors = self.vector.reshape(1, -1)[0]
        paths = self.path.reshape(1, -1)[0]
        result = np.concatenate((points, vectors, paths, self.number, self.wavelength, self.incidence, self.intensity))
        result = result[~np.isnan(result)]
        result.flags.writeable = False
        result_hash = hash(result.data.tobytes())
        return result_hash

    def __repr__(self) -> str:
        intro = f"RayList with {self.N} rays\n"
        average = mgeo.Vector(np.mean(self.vector)).normalized()
        vectors = self.vector.normalized()
        angles = np.arccos(np.clip(np.dot(vectors, average), -1.0, 1.0))
        avg_string = "On average, the rays are oriented along the vector:\n" + str(average) + "\n"
        NA = f"Numerical aperture: {np.max(np.sin(angles)):.3f}\n"
        return intro +avg_string+ NA