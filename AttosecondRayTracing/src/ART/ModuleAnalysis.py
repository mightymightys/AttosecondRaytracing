"""
Adds various analysis methods.

Created in July 2024

@author: André Kalouguine + Stefan Haessler + Anthony Guillaume
"""
# %% Module imports
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleProcessing as mp
import ART.ModulePlottingMethods as mpm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math

LightSpeed = 299792458000

# %% Analysis methods
def getETransmission(RayListIn, RayListOut) -> float:
    """
    Calculates the energy transmission from RayListIn to RayListOut in percent by summing up the
    intensity-property of the individual rays.

    Parameters
    ----------
        RayListIn : list(Ray)
            List of incoming rays.
        
        RayListOut : list(Ray)
            List of outgoing rays.

    Returns
    -------
        ETransmission : float
    """
    ETransmission = 100 * sum(Ray.intensity for Ray in RayListOut) / sum(Ray.intensity for Ray in RayListIn)
    return ETransmission


def GetResultSummary(Detector, RayListAnalysed, verbose=False):
    """
    Calculate and return FocalSpotSize-standard-deviation and Duration-standard-deviation
    for the given Detector and RayList.
    If verbose, then also print a summary of the results for the given Detector.

    Parameters
    ----------
        Detector : Detector
            An object of the ModuleDetector.Detector-class.

        RayListAnalysed : list(Ray)
            List of objects of the ModuleOpticalRay.Ray-class.

        verbose : bool
            Whether to print a result summary.

    Returns
    -------
        FocalSpotSizeSD : float

        DurationSD : float
    """
    DetectorPointList2DCentre = Detector.get_PointList2DCentre(RayListAnalysed)
    FocalSpotSizeSD = mp.StandardDeviation(DetectorPointList2DCentre)
    DelayList = Detector.get_Delays(RayListAnalysed)
    DurationSD = mp.StandardDeviation(DelayList)

    if verbose:
        FocalSpotSize = mgeo.DiameterPointList(DetectorPointList2DCentre)
        summarystring = (
            "At the detector distance of "
            + "{:.3f}".format(Detector.get_distance())
            + " mm we get:\n"
            + "Spatial std : "
            + "{:.3f}".format(FocalSpotSizeSD * 1e3)
            + " \u03BCm and min-max: "
            + "{:.3f}".format(FocalSpotSize * 1e3)
            + " \u03BCm\n"
            + "Temporal std : "
            + "{:.3e}".format(DurationSD)
            + " fs and min-max : "
            + "{:.3e}".format(max(DelayList) - min(DelayList))
            + " fs"
        )

        print(summarystring)

    return FocalSpotSizeSD, DurationSD


def ReturnNumericalAperture(RayList: list[mray.Ray], RefractiveIndex: float = 1) -> float:
    r"""
    Returns the numerical aperture associated with the supplied ray-bundle 'Raylist'.
    This is $n\sin\theta$, where $\theta$ is the maximum angle between any of the rays and the central ray,
    and $n$ is the refractive index of the propagation medium.

    Parameters
    ----------
        RayList : list of Ray-object
            The ray-bundle of which to determine the NA.

        RefractiveIndex : float, optional
            Refractive index of the propagation medium, defaults to =1.

    Returns
    -------
        NA : float
    """
    CentralRay = mp.FindCentralRay(RayList)
    if CentralRay is None:
        CentralVector = np.array([0, 0, 0])
        for k in RayList:
            CentralVector = CentralVector + k.vector
        CentralVector = CentralVector / len(RayList)
    else:
        CentralVector = CentralRay.vector
    ListAngleAperture = []
    for k in RayList:
        ListAngleAperture.append(mgeo.AngleBetweenTwoVectors(CentralVector, k.vector))

    return np.sin(np.amax(ListAngleAperture)) * RefractiveIndex


def ReturnAiryRadius(Wavelength: float, NumericalAperture: float) -> float:
    r"""
    Returns the radius of the Airy disk: $r = 1.22 \frac\{\lambda\}\{NA\}$,
    i.e. the diffraction-limited radius of the focal spot corresponding to a given
    numerical aperture $NA$ and a light wavelength $\lambda$.

    For very small $NA<10^\{-3\}$, diffraction effects becomes negligible and the Airy Radius becomes meaningless,
    so in that case, a radius of 0 is returned.

    Parameters
    ----------
        Wavelength : float
            Light wavelength in mm.

        NumericalAperture : float

    Returns
    -------
        AiryRadius : float
    """
    if NumericalAperture > 1e-3 and Wavelength is not None:
        return 1.22 * 0.5 * Wavelength / NumericalAperture
    else:
        return 0  # for very small numerical apertures, diffraction effects becomes negligible and the Airy Radius becomes meaningless


def get_planewavefocus(OpticalChain, Detector, Index, size=None, Nrays=1000, resolution=100):
    """
    This function calculates the approximate polychromatic focal spot of a set of rays.
    To do so, it positions itself in the detector plane and samples a square area.
    If the size of the area is not given, it will use twice the Airy radius of the system with the largest wavbelength.
    Otherwise it uses the size.
    The resolution of the sampling is given by the resolution parameter. So it samples on a grid resolution x resolution.

    To calculate the intensity, it takes Nrays out of the raylist (to subsample if needed).
    It assimilates every ray to a plane wave, so it calculates a k-vector for each ray (taking into account the wavelength of the ray).
    It calculates the phase of each ray from the delay and from the intersection position with the detector.
    The delay from the non-central position is simply sin(alpha)*distance/c, where alpha is the angle between the ray and the normal to the detector
    and distance is the distance between the intersection point and the central point of the detector. 
    It then calculates the intensity at each point of the grid by summing the intensity of each plane wave.

    As long as there are not too many rays, this method is faster than doing an FFT.
    It's also naturally polychromatic.
    """
    RayList = OpticalChain.get_output_rays()[Index]
    if size is None:
        Wavelengths = [Ray.wavelength for Ray in RayList]
        Wavelength = max(Wavelengths)
        NumericalAperture = ReturnNumericalAperture(RayList)
        size = 3 * ReturnAiryRadius(Wavelength, NumericalAperture)
    X = np.linspace(-size / 2, size / 2, resolution)
    Y = np.linspace(-size / 2, size / 2, resolution)
    X, Y = np.meshgrid(X, Y)
    # We pick Nrays. if there are less than Nrays, we take all of them.
    PickedRaysGlobal = np.random.choice(RayList, min(Nrays, len(RayList)), replace=False)
    # We calculate the k-vector of each ray, taking into account the wavelength of the ray.
    # The units should be in mm^-1
    PickedRays = [Ray.to_basis(*Detector.basis) for Ray in PickedRaysGlobal]
    # The rays are now in the reference plane of the detector whose normal is [0,0,1]
    wavelengths = np.array([Ray.wavelength for Ray in PickedRays])
    vectors = np.array([Ray.vector for Ray in PickedRays])
    frequencies = np.array([LightSpeed / Ray.wavelength for Ray in PickedRays]) # mm/s / mm = 1/s
    k_vectors = 2*np.pi * vectors / wavelengths[:, np.newaxis]
    angles = np.arccos(np.clip(np.dot(vectors, np.array([0, 0, 1])), -1.0, 1.0))
    # We calculate the intersection of the rays with the detector plane
    Intersections = mgeo.IntersectionRayListZPlane(PickedRays)[0]
    distances = (Intersections - mgeo.Origin[:2]).norm
    # We calculate the phase of each ray
    PathDelays = np.array(Detector.get_Delays(PickedRaysGlobal))
    PositionDelays = np.sum(np.array(Intersections._add_dimension()-mgeo.Origin)*vectors, axis=1)  / LightSpeed * 1e15
    Delays = (PathDelays+PositionDelays) /  1e15
    # We have the delays in fs, we can now calculate the phase, taking into account the wavelength of each ray
    Phases = np.array([2 * np.pi * frequencies[i] * Delays[i] for i in range(len(Delays))])
    # We also need the intensity of each ray
    RayIntensities = np.array([Ray.intensity for Ray in PickedRays])
    # We can now calculate the intensity at each point of the grid
    Intensity = np.zeros((resolution, resolution))
    for i in range(len(PickedRays)):
        Intensity += RayIntensities[i] * np.cos(k_vectors[i][0] * X + k_vectors[i][1] * Y + Phases[i])

    return X, Y, Intensity


def get_diffractionfocus(OpticalChain, Detector, Index, size=None, Nrays=1000, resolution=100):
    """
    This function calculates the approximate polychromatic focal spot of a set of rays.
    To do so, it positions itself in the detector plane and samples a square area.
    If the size of the area is not given, it will use twice the Airy radius of the system with the largest wavbelength.
    Otherwise it uses the size.
    The resolution of the sampling is given by the resolution parameter. So it samples on a grid resolution x resolution.

    To calculate the intensity, it takes Nrays out of the raylist (to subsample if needed).
    It assimilates every ray to a plane wave, so it calculates a k-vector for each ray (taking into account the wavelength of the ray).
    It considers that all the rays are intersecting the detector in the middle and doesn't take into account their phase.

    So it returns the best case scenario for a diffraction limited focus with that numerical aperture
    """
    RayList = OpticalChain.get_output_rays()[Index]
    if size is None:
        Wavelengths = [Ray.wavelength for Ray in RayList]
        Wavelength = max(Wavelengths)
        NumericalAperture = ReturnNumericalAperture(RayList)
        size = 3 * ReturnAiryRadius(Wavelength, NumericalAperture)
    X = np.linspace(-size / 2, size / 2, resolution)
    Y = np.linspace(-size / 2, size / 2, resolution)
    X, Y = np.meshgrid(X, Y)
    # We pick Nrays. if there are less than Nrays, we take all of them.
    PickedRaysGlobal = np.random.choice(RayList, min(Nrays, len(RayList)), replace=False)
    # We calculate the k-vector of each ray, taking into account the wavelength of the ray.
    # The units should be in mm^-1
    PickedRays = [Ray.to_basis(*Detector.basis) for Ray in PickedRaysGlobal]
    # The rays are now in the reference plane of the detector whose normal is [0,0,1]
    wavelengths = np.array([Ray.wavelength for Ray in PickedRays])
    vectors = np.array([Ray.vector for Ray in PickedRays])
    frequencies = np.array([LightSpeed / Ray.wavelength for Ray in PickedRays]) # mm/s / mm = 1/s
    k_vectors = 2*np.pi * vectors / wavelengths[:, np.newaxis]
    angles = np.arccos(np.clip(np.dot(vectors, np.array([0, 0, 1])), -1.0, 1.0))
    # We calculate the intersection of the rays with the detector plane
    Intersections = mgeo.IntersectionRayListZPlane(PickedRays)[0]
    distances = (Intersections - mgeo.Origin[:2]).norm
    # We also need the intensity of each ray
    RayIntensities = np.array([Ray.intensity for Ray in PickedRays])
    # We can now calculate the intensity at each point of the grid
    Intensity = np.zeros((resolution, resolution))
    for i in range(len(PickedRays)):
        Intensity += RayIntensities[i] * np.cos(k_vectors[i][0] * X + k_vectors[i][1] * Y)

    return X, Y, Intensity

# %% Asphericity analysis
# This code is to do asphericity analysis of various surfaces
# It calculates the closes sphere to the surface of an optical element
# Then there are two functions. One that simply gives an asphericity value
# and another one that actually plots the distance the the closest sphere 
# in much the same way as we plot the MirrorProjection

def BestFitSphere(X,Y,Z):
    """
    This function calculates the best sphere to fit a set of points.
    It uses the least square method to find the center and the radius of the sphere.
    Cite https://jekel.me/2015/Least-Squares-Sphere-Fit/
    """
    A = np.zeros((len(X),4))
    A[:,0] = X*2
    A[:,1] = Y*2
    A[:,2] = Z*2
    A[:,3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(X),1))
    f[:,0] = X**2 + Y**2 + Z**2
    C, residules, rank, singval = np.linalg.lstsq(A,f)

    #   solve for the radius
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = math.sqrt(t)

    return mgeo.Point(C[:3].flatten()),radius


def get_closest_sphere(Mirror, Npoints=1000):
    """
    This function calculates the closest sphere to the surface of a mirror.
    It does so by sampling the surface of the mirror at Npoints points.
    It then calculates the closest sphere to these points.
    It returns the radius of the sphere and the center of the sphere.
    """
    Points = mpm.sample_support(Mirror.support, Npoints=1000)
    Points += Mirror.r0[:2]
    Z = Mirror._zfunc(Points)
    Points = mgeo.PointArray([Points[:, 0], Points[:, 1], Z]).T
    spX, spY, spZ = Points[:, 0], Points[:, 1], Points[:, 2]
    Center, Radius = BestFitSphere(spX, spY, spZ)
    return Center, Radius

def get_asphericity(Mirror, Npoints=1000):
    """
    This function calculates the maximum distance of the mirror surface to the closest sphere. 
    """
    center, radius = get_closest_sphere(Mirror, Npoints)
    Points = mpm.sample_support(Mirror.support, Npoints=1000)
    Points += Mirror.r0[:2]
    Z = Mirror._zfunc(Points)
    Points = mgeo.PointArray([Points[:, 0], Points[:, 1], Z]).T
    Points_centered = Points - center
    Distance = np.linalg.norm(Points_centered, axis=1) - radius
    Distance*=1e3 # To convert to µm
    return np.ptp(Distance)

