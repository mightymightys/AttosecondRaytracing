"""
Provides classes for different deformations of optical surfaces.


Created on Fri Mar 10 2023

@author: Semptum + Stefan Haessler 
"""
# %% Modules

import numpy as np
from abc import ABC, abstractmethod
import scipy
import math
from ART.recursive_zernike_generator import zernike_gradient

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
        rect = Support._CircumRect()
        X = np.linspace(-rect[0], rect[0], num=self.deformation.shape[0])
        Y = np.linspace(-rect[1], rect[1], num=self.deformation.shape[1])
        self.DerivX, self.DerivY = np.gradient(self.deformation, rect/self.deformation.shape)
        DerivInterpX = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.DerivX), method="linear")
        DerivInterpY = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.DerivY), method="linear")
        SurfInterp = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(self.deformation), method="linear")

        self.rms = np.std(Map)
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
    def __init__(self, Support, RMS, slope=-2, smallest=0.01, biggest=10):
        # The sizes are the wavelength in mm
        rect = Support._CircumRect()
        k_max = np.pi / smallest
        k_min = np.pi / biggest
        ResX = int(2 ** math.ceil(math.log2(k_max * rect[0] / 2)))
        ResY = int(2 ** math.ceil(math.log2(k_max * rect[1] / 2)))
        ResX = int(round(k_max * rect[0] / 2))
        ResY = int(round(k_max * rect[1] / 2))

        # TODO Ensure that we don't get astigmatic pixels because of log/ceil
        if slope == -2:
            integ = np.log(k_max / k_min)
        else:
            a = slope + 2
            integ = (k_max**a - k_min**a) / a
        # kX = np.linspace(0, k_max, num=ResX)
        # kY = np.linspace(0, k_max, num=ResY)
        # kXX, kYY = np.meshgrid(kX, kY)
        kXX, kYY = np.meshgrid(
            np.linspace(0, k_max, num=ResX, dtype='float32'),
            np.linspace(0, k_max, num=ResY, dtype='float32'),
            sparse=True)

        ## 1stversion
        FFT = np.sqrt(kXX**2 + kYY**2, dtype='complex64')  # Matrix of k norms
        mask = (FFT > k_min) & (FFT < k_max)
        FFT[np.logical_not(mask)] = 0.0
        FFT[mask] = np.sqrt(FFT[mask]**slope * (2 / np.pi * RMS**2 / integ))
        del mask

        # ## 2nd version
        # D = np.sqrt(kXX**2 + kYY**2)  # Matrix of k norms
        # mask = (D > k_min) & (D < k_max)
        # #b = np.log(2 / np.pi * RMS**2 / integ)
        # FFT = np.zeros_like(D)
        # #FFT[mask] = D[mask] ** slope * np.exp(b) #why do we first calc. b as a log and then use its exp()? Just leave both away? 
        # FFT[mask] = D[mask] **slope * (2 / np.pi * RMS**2 / integ)
        # FFT = np.sqrt(FFT)

        def tiled(x):
            return np.block([
                [np.transpose(np.rot90(x)), x],
                [np.rot90(x, k=2), np.transpose(np.rot90(x, k=-1))]
            ])
        phase = np.random.uniform(0, 2 * np.pi, size=FFT.shape).astype('float32')
        FFT = FFT * np.exp(1j * phase)
        del phase
        FFT_tiled = tiled(FFT)
        del FFT
        gridsize = FFT_tiled.shape

        deformation = np.fft.irfft2(np.fft.fftshift(FFT_tiled), gridsize)

        kXX2, kYY2 = np.meshgrid(np.linspace(-k_max, k_max, num=ResX * 2,  dtype='float32'), np.linspace(-k_max, k_max, num=ResY * 2,  dtype='float32'), sparse=True)
        #DerivFFTX_tiled = FFT_tiled * 1j * kXX2
        #DerivFFTY_tiled = FFT_tiled * 1j * kYY2 
        DerivX = np.fft.irfft2(np.fft.fftshift(FFT_tiled * 1j * kXX2), gridsize)
        DerivY = np.fft.irfft2(np.fft.fftshift(FFT_tiled * 1j * kYY2), gridsize)
        del FFT_tiled

        X = np.linspace(-rect[0], rect[0], num=ResX * 2)  # TODO why
        Y = np.linspace(-rect[1], rect[1], num=ResY * 2)

        DerivInterpX = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(DerivX), method="linear")
        DerivInterpY = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(DerivY), method="linear")
        SurfInterp = scipy.interpolate.RegularGridInterpolator((X, Y), np.transpose(deformation), method="linear")

        self.DerivX = DerivX
        self.DerivY = DerivY
        self.rms = np.sqrt(np.mean(deformation**2))
        self.deformation = deformation
        self.DerivInterp = lambda x: np.array([DerivInterpX(x[:2]), DerivInterpY(x[:2])])
        self.SurfInterp = lambda x: SurfInterp(x[:2])

        print(self.rms)

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


def normal_add(N1, N2):
    """
    Simple function that takes in two normal vectors of a deformation and calculates the total normal vector if the two deformations were individually applied.
    """
    normal1 = N1 / np.linalg.norm(N1)
    normal2 = N2 / np.linalg.norm(N2)
    grad1X = -normal1[0] / normal1[2]
    grad1Y = -normal1[1] / normal1[2]
    grad2X = -normal2[0] / normal2[2]
    grad2Y = -normal2[1] / normal2[2]
    gradX = grad1X + grad2X
    gradY = grad1Y + grad2Y
    return np.array([-gradX, -gradY, 1])
