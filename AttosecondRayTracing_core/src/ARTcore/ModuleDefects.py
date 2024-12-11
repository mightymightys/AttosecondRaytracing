"""
Provides classes for different deformations of optical surfaces.


Created on Fri Mar 10 2023

@author: Semptum + Stefan Haessler 
"""
# %% Modules

import numpy as np
#import cupy as cp
from abc import ABC, abstractmethod
import scipy
import math
from ARTcore.ModuleZernike import zernike_gradient
import logging

logger = logging.getLogger(__name__)

# %% Abstract class


class Defect(ABC):
    """
    Abstract class representing a defect on the optical surface.
    """

    @abstractmethod
    def RMS(self):
        pass

    @abstractmethod
    def PV(self):
        pass

class MeasuredMap(Defect):
    """
    Class describing a defect map measured experimentally (for instance using a Hartmann wavefront sensor). Note that the provided image should cover the entire support. 
    """
    def __init__(self, Support, Map):
        self.deformation = Map
        self. Support = Support
        rect = Support._CircumRect() # TODO fix this so it just uses a rectangle with known size
        X = np.linspace(-rect[0], rect[0], num=self.deformation.shape[0])
        Y = np.linspace(-rect[1], rect[1], num=self.deformation.shape[1])
        self.DerivX, self.DerivY = np.gradient(self.deformation, rect/self.deformation.shape)
        DerivInterpX = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.DerivX), method="linear")
        DerivInterpY = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.DerivY), method="linear")
        SurfInterp = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.deformation), method="linear")

        self.rms = np.std(self.deformation)
        self.DerivInterp = lambda x: np.array([DerivInterpX(x[:2]), DerivInterpY(x[:2])])
        self.SurfInterp = lambda x: SurfInterp(x[:2])
    def get_normal(self, Point):
        Grad = self.DerivInterp(Point)
        dX, dY =  Grad.flatten()
        norm = np.linalg.norm([dX, dY, 1])
        dX /= norm
        dY /= norm
        return np.array([dX, dY, np.sqrt(1 - dX**2 - dY**2)])

    def get_offset(self, Point):
        return self.SurfInterp(Point)

    def RMS(self):
        return self.rms

    def PV(self):
        pass

class Fourrier(Defect):
    """
    
    """
    def __init__(self, Support, RMS, slope=-2, smallest=0.1, biggest=None):
        # The sizes are the wavelength in mm
        rect = Support._CircumRect()  # TODO fix this so it just uses a rectangle with known size
        if biggest is None:
            biggest = np.max(rect)
            
        k_max = 2 / smallest
        k_min = 2 / biggest
        ResX = int(round(k_max * rect[0] / 2))+1
        ResY = int(round(k_max * rect[1]))

        kXX, kYY = np.meshgrid(
            np.linspace(0, k_max, num=ResX, dtype='float32', endpoint=False),
            np.linspace(-k_max, k_max, num=ResY, dtype='float32', endpoint=False),
            sparse=True)

        maskedFFT = np.ma.masked_outside(np.sqrt(kXX**2 + kYY**2), k_min, k_max)

        FFT = maskedFFT**slope * np.exp(
            1j *np.random.uniform(0, 2 * np.pi, size=maskedFFT.shape).astype('float32')
            )
        FFT = FFT.data*(1-FFT.mask)

        deformation = np.fft.irfft2(np.fft.ifftshift(FFT, axes=0))
        RMS_factor = RMS/np.std(deformation)
        deformation *= RMS_factor

        DerivX = np.fft.irfft2(np.fft.ifftshift(FFT * 1j * kXX * RMS_factor, axes=0))*np.pi/2
        kY = np.concatenate((kYY[kYY.shape[0]//2:],kYY[:kYY.shape[0]//2]))
        DerivY = np.fft.irfft2(np.fft.ifftshift(FFT * 1j * RMS_factor, axes=0)*kY)*np.pi/2
        #del FFT

        X = np.linspace(-rect[0]/2, rect[0]/2, num=(ResX-1)*2) # Because iRfft
        Y = np.linspace(-rect[1]/2, rect[1]/2, num=ResY)
        
        DerivInterpX = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(DerivX), method="linear")
        DerivInterpY = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(DerivY), method="linear")
        SurfInterp = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(deformation), method="linear")

        self.DerivX = DerivX
        self.DerivY = DerivY
        self.rms = np.std(deformation)
        self.deformation = deformation
        self.DerivInterp = lambda x: np.array([DerivInterpX(x[:2]), DerivInterpY(x[:2])])
        self.SurfInterp = lambda x: SurfInterp(x[:2])

    def get_normal(self, Point):
        """
        Calculates and returns the surface normal at the given Point.
        It uses the derivative interpolators to compute the partial derivatives of the surface
        deformation and returns a normalized vector representing the surface normal.
        """
        dX, dY = self.DerivInterp(Point)
        norm = np.linalg.norm([dX, dY, 1])
        dX /= norm
        dY /= norm
        return np.array([dX, dY, np.sqrt(1 - dX**2 - dY**2)])

    def get_offset(self, Point):
        """
        Calculates and returns the surface offset at the given Point.
        It uses the surface deformation interpolator to determine the displacement
        at the specified point on the surface.
        """
        return self.SurfInterp(Point)

    def RMS(self):
        """
        Returns the root mean square (RMS) value of the surface deformation (self.rms).
        """
        return self.rms

    def PV(self):
        pass


class Zernike(Defect):
    def __init__(self, Support, coefficients):
        self.coefficients = coefficients
        self.max_order = np.max([k[0] for k in coefficients])
        self.support = Support
        self.R = Support._CircumCirc()

    def get_normal(self, Point):
        x, y, z = Point / self.R
        val, derivX, derivY = zernike_gradient([x], [y], self.max_order)
        dX = 0
        dY = 0
        for k in self.coefficients:
            dX += self.coefficients[k] * derivX[k][0][1][0]
            dY += self.coefficients[k] * derivY[k][0][1][0]
        dX /= self.R
        dY /= self.R
        return np.array([-dX, -dY, 1])

    def get_offset(self, Point):
        x, y, z = Point / self.R
        val, derivX, derivY = zernike_gradient([x], [y], self.max_order)
        Z = 0
        for k in self.coefficients:
            Z += self.coefficients[k] * val[k][0][1][0]
        return Z

    def RMS(self):
        return np.sqrt(np.sum([i**2 for i in self.coefficients.values()]))

    def PV(self):
        pass

