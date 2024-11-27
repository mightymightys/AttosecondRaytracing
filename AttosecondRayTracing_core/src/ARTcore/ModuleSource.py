"""
Provides a definition of sources of light.
A source can be either a composite source (a sum of sources) or a simple source.
A simple source is defined by 4 parameters:
- Spectrum
- Angular power distribution
- Ray origins distribution
- Ray directions distribution
A reason to use composite sources is for when the power and ray distributions depend on the wavelength.
Another source option is the source containing just a list of rays that can be prepared manually.
Finally another source option is one used for alignment, containing a single ray.


Created in 2019

@author: Anthony Guillaume and Stefan Haessler
"""

# %% Modules
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin
import ARTcore.ModuleProcessing as mp

from ARTcore.DepGraphDefinitions import UniformSpectraCalculator

import numpy as np
from abc import ABC, abstractmethod
import logging

logger = logging.getLogger(__name__)

# %% Abstract base classes for sources
class Source(ABC):
    """
    Abstract base class for sources.
    Fundamentally a source just has to be able to return a list of rays.
    To do so, we make it callable.
    """
    @abstractmethod
    def __call__(self,N):
        pass

class PowerDistribution(ABC):
    """
    Abstract base class for angular power distributions.
    """
    @abstractmethod
    def __call__(self, Origins, Directions):
        """
        Return the power of the rays coming 
        from the origin points Origins to the directions Directions.
        """
        pass
class RayOriginsDistribution(ABC):
    """
    Abstract base class for ray origins distributions.
    """
    @abstractmethod
    def __call__(self, N):
        """
        Return the origins of N rays.
        """
        pass

class RayDirectionsDistribution(ABC):
    """
    Abstract base class for ray directions distributions.
    """
    @abstractmethod
    def __call__(self, N):
        """
        Return the directions of N rays.
        """
        pass

class Spectrum(ABC):
    """
    Abstract base class for spectra.
    """
    @abstractmethod
    def __call__(self, N):
        """
        Return the wavelengths of N rays.
        """
        pass

# %% Specific power distributions
class SpatialGaussianPowerDistribution(PowerDistribution):
    """
    Spatial Gaussian power distribution.
    """
    def __init__(self, Power, W0):
        self.Power = Power
        self.W0 = W0

    def __call__(self, Origins, Directions):
        """
        Return the power of the rays coming 
        from the origin points Origins to the directions Directions.
        """
        return self.Power * np.exp(-2 * np.linalg.norm(Origins, axis=1) ** 2 / self.W0 ** 2)
    
class AngularGaussianPowerDistribution(PowerDistribution):
    """
    Angular Gaussian power distribution.
    """
    def __init__(self, Power, Divergence):
        self.Power = Power
        self.Divergence = Divergence

    def __call__(self, Origins, Directions):
        """
        Return the power of the rays coming 
        from the origin points Origins to the directions Directions.
        """
        return self.Power * np.exp(-2 * np.arccos(np.dot(Directions, [0, 0, 1])) ** 2 / self.Divergence ** 2)

class GaussianPowerDistribution(PowerDistribution):
    """
    Gaussian power distribution.
    """
    def __init__(self, Power, W0, Divergence):
        self.Power = Power
        self.W0 = W0
        self.Divergence = Divergence

    def __call__(self, Origins, Directions):
        """
        Return the power of the rays coming 
        from the origin points Origins to the directions Directions.
        """
        return self.Power * np.exp(-2 * (Origins-mgeo.Origin).norm ** 2 / self.W0 ** 2) * np.exp(-2 * np.array(np.arccos(np.dot(Directions, mgeo.Vector([1, 0, 0])))) ** 2 / self.Divergence ** 2)

class UniformPowerDistribution(PowerDistribution):
    """
    Uniform power distribution.
    """
    def __init__(self, Power):
        self.Power = Power

    def __call__(self, Origins, Directions):
        """
        Return the power of the rays coming 
        from the origin points Origins to the directions Directions.
        """
        return self.Power * np.ones(len(Origins))

# %% Specific ray origins distributions
class PointRayOriginsDistribution(RayOriginsDistribution):
    """
    Point ray origins distribution.
    """
    def __init__(self, Origin):
        self.Origin = Origin

    def __call__(self, N):
        """
        Return the origins of N rays.
        """
        return mgeo.PointArray([self.Origin for i in range(N)])
    
class DiskRayOriginsDistribution(RayOriginsDistribution):
    """
    Disk ray origins distribution. Uses the Vogel spiral to initialize the rays.
    """
    def __init__(self, Origin, Radius):
        self.Origin = Origin
        self.Radius = Radius

    def __call__(self, N):
        """
        Return the origins of N rays.
        """
        MatrixXY = mgeo.SpiralVogel(N, self.Radius)
        return mgeo.PointArray([self.Origin + mgeo.Vector([MatrixXY[i, 0], MatrixXY[i, 1], 0]) for i in range(N)])
    
# %% Specific ray directions distributions
class UniformRayDirectionsDistribution(RayDirectionsDistribution):
    """
    Uniform ray directions distribution.
    """
    def __init__(self, Direction):
        self.Direction = Direction

    def __call__(self, N):
        """
        Return the directions of N rays.
        """
        return mgeo.VectorArray([self.Direction for i in range(N)])

class ConeRayDirectionsDistribution(RayDirectionsDistribution):
    """
    Cone ray directions distribution. Uses the Vogel spiral to initialize the rays.
    """
    def __init__(self, Direction, Angle):
        self.Direction = Direction
        self.Angle = Angle

    def __call__(self, N):
        """
        Return the directions of N rays.
        """
        Height = 1
        Radius = Height * np.tan(self.Angle)
        MatrixXY = mgeo.SpiralVogel(N, Radius)
        q = mgeo.QRotationVector2Vector(mgeo.Vector([0, 0, 1]), self.Direction)
        return mgeo.VectorArray([[MatrixXY[i, 0], MatrixXY[i, 1], Height] for i in range(N)]).rotate(q).normalized()


# %% Specific spectra
class SingleWavelengthSpectrum(Spectrum):
    """
    Single wavelength spectrum.
    """
    def __init__(self, Wavelength):
        self.Wavelength = Wavelength

    def __call__(self, N):
        """
        Return the wavelengths of N rays.
        """
        return np.ones(N) * self.Wavelength

class UniformSpectrum(Spectrum):
    """
    Uniform spectrum.
    Can be specified as a r
    """
    def __init__(self, eVMax = None, eVMin = None, eVCentral = None, 
                 lambdaMax = None, lambdaMin = None, lambdaCentral = None, 
                 eVWidth = None, lambdaWidth = None):
        # Using the DepSolver, we calculate the minimum and maximum wavelengths
        values, steps = UniformSpectraCalculator().calculate_values(
            lambda_min = lambdaMin, 
            lambda_max = lambdaMax, 
            lambda_center = lambdaCentral, 
            lambda_width = lambdaWidth,
            eV_min = eVMin,
            eV_max = eVMax,
            eV_center = eVCentral,
            eV_width = eVWidth
        )
        self.lambdaMin = values['lambda_min']
        self.lambdaMax = values['lambda_max']
    def __call__(self, N):
        """
        Return the wavelengths of N rays.
        """
        return np.linspace(self.lambdaMin, self.lambdaMax, N)


# %% Simple sources
class SimpleSource(Source):
    """
    A simple source is defined by 4 parameters:
    - Spectrum
    - Power distribution
    - Ray origins distribution
    - Ray directions distribution
    """

    def __init__(self, Spectrum, PowerDistribution, RayOriginsDistribution, RayDirectionsDistribution):
        self.Spectrum = Spectrum
        self.PowerDistribution = PowerDistribution
        self.RayOriginsDistribution = RayOriginsDistribution
        self.RayDirectionsDistribution = RayDirectionsDistribution

    def __call__(self, N):
        """
        Return a list of N rays from the simple source.
        """
        Wavelengths = self.Spectrum(N)
        Origins = self.RayOriginsDistribution(N)
        Directions = self.RayDirectionsDistribution(N)
        Powers = self.PowerDistribution(Origins, Directions)
        RayList = []
        for i in range(N):
            RayList.append(mray.Ray(Origins[i], Directions[i], wavelength=Wavelengths[i], number=i, intensity=Powers[i]))
        return mray.RayList.from_list(RayList)


class ListSource(Source):
    """
    A source containing just a list of rays that can be prepared manually.
    """
    def __init__(self, Rays):
        self.Rays = Rays

    def __call__(self, N):
        """
        Return the list of rays.
        """
        if N < len(self.Rays):
            return self.Rays[:N]
        elif N>len(self.Rays):
            logger.warning("Requested number of rays is greater than the number of rays in the source. Returning the whole source.")
            return self.Rays    
        return mray.RayList.from_list(self.Rays)