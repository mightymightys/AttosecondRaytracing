"""
Provides functions that create lists of Ray-objects that can serve as light sources for the start of the ray tracing calculation.

These are called by 'ARTmain' according to the values in the *SourceProperties*-dictionary filled in the CONFIG-scripts..



Created in 2019

@author: Anthony Guillaume and Stefan Haessler
"""

# %% Modules

import numpy as np
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleProcessing as mp

# %%


def _Cone(Angle: float, NbRays: int, Wavelength=None):
    """
    Return a list of rays filling a cone according to Vogel's spiral.

    Parameters
    ----------
        Angle : float
            Cone apex *half*-angle in radians.

        NbRays : int
            Number of rays in the list.

        Wavelength : float, optional
            Wavelength of the light, set as attribute to all rays.
    """
    Height = 1
    Radius = Height * np.tan(Angle)

    MatrixXY = mgeo.SpiralVogel(NbRays, Radius)

    RayList = []

    for k in range(NbRays):
        x = MatrixXY[k, 0]
        y = MatrixXY[k, 1]
        RayList.append(mray.Ray(np.array([0, 0, 0]), np.array([x, y, Height]), Number=k, Wavelength=Wavelength))

    return RayList


# %%
def PointSource(S: np.ndarray, Axis: np.ndarray, Divergence: float, NbRays: int, Wavelength=None):
    """
    Return a list of rays simulating rays from a point source.

    ![Illustration of a point source.](pointsource.svg)

    Parameters
    ----------
        S : np.ndarray
            Coordinate vector of the source point (in mm).

        Axis : np.ndarray
            Axis vector of the cone (main direction of the rays)

        Divergence : float
            Cone apex *half*-angle in *radians* (similar to gaussian beam).

        NbRays : int
            Number of rays in the list.

        Wavelength : float, optional
            Wavelength of the light, set as attribute to all rays.

    """
    RayList = _Cone(Divergence, NbRays, Wavelength=Wavelength)
    RayList = mgeo.RotationRayList(RayList, np.array([0, 0, 1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, S)
    return RayList


# %%
def ExtendedSource(S: np.ndarray, Axis: np.ndarray, Diameter: float, Divergence: float, NbRays: int, Wavelength=None):
    """
    Return a list of rays simulating rays from an extended array of point sources, distributed over a disk with Diameter.

    Parameters
    ----------
        S : np.ndarray
            Coordinate vector of the source point (in mm).

        Axis : np.ndarray
            Axis vector of the cone (main direction of the rays)

        Diameter : float
            Diameter of the extended source

        Divergence : float
            Cone apex *half*-angle in *radians* (similar to gaussian beam).

        NbRays : int
            Number of rays in the list.

        Wavelength : float, optional
            Wavelength of the light, set as attribute to all rays.

    """
    NbPointSources = max(5, int(250 * Diameter))
    NbPointSources = min(NbPointSources, 100)

    MatrixXY = mgeo.SpiralVogel(NbPointSources, Diameter / 2)  # the positions of the point sources

    NbRaysPerPointSource = max(100, int(NbRays / NbPointSources))
    RayList = []
    PointSourceRayList = _Cone(Divergence, NbRaysPerPointSource, Wavelength=Wavelength)
    for k in range(NbPointSources):
        ShiftedPointSourceRayList = mgeo.TranslationRayList(PointSourceRayList, [MatrixXY[k, 0], MatrixXY[k, 1], 0])
        for l in range(NbRaysPerPointSource):
            ShiftedPointSourceRayList[l]._number += (
                k * NbRaysPerPointSource
            )  # we allow ourselves exceptionally to modify the "internal" _number attribute
        RayList.extend(ShiftedPointSourceRayList)

    RayList = mgeo.RotationRayList(RayList, np.array([0, 0, 1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, S)
    return RayList


# %%
def PlaneWaveDisk(Centre: np.ndarray, Axis: np.ndarray, Radius: float, NbRays: int, Wavelength=None):
    """
    Return a list of rays, all parallel and distributed over a round cross section according to a Vogel-spiral (see illustration). Simulating a 'round' collimated beam.

    ![Illustration of a Vogel spiral.](doc_VogelSpiral.png)

    Parameters
    ----------
        Centre : np.ndarray
            Coordinate vector of the center of the source disk (in mm).

        Axis : np.ndarray
            Axis vector of the the plane wave (vector of all rays)

        Radius : float
            Radius of the source ray-bundle.

        NbRays : int
            Number of rays in the list.

        Wavelength : float, optional
            Wavelength of the light, set as attribute to all rays.

    """
    MatrixXY = mgeo.SpiralVogel(NbRays, Radius)
    RayList = []
    num = 0
    for k in range(NbRays - 1):
        x = MatrixXY[k, 0]
        y = MatrixXY[k, 1]
        RayList.append(mray.Ray(np.array([x, y, 0]), np.array([0, 0, 1]), Number=num, Wavelength=Wavelength))
        num = num + 1
    RayList = mgeo.RotationRayList(RayList, np.array([0, 0, 1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, Centre)
    return RayList


# %%
def PlaneWaveSquare(Centre: np.ndarray, Axis: np.ndarray, SideLength: float, NbRays: int, Wavelength=None):
    """
    Return a list of rays, all parallel and distributed over a square cross section awith SideLength. Simulating a 'square' collimated beam.

    Parameters
    ----------
        Centre : np.ndarray
            Coordinate vector of the center of the source disk (in mm).

        Axis : np.ndarray
            Axis vector of the the plane wave (vector of all rays)

        SideLength : float
            Side length of the square source ray-bundle.

        NbRays : int
            Number of rays in the list.

        Wavelength : float, optional
            Wavelength of the light, set as attribute to all rays.

    """
    RayList = []
    x = np.linspace(-SideLength / 2, SideLength / 2, int(np.sqrt(NbRays)))
    y = np.linspace(-SideLength / 2, SideLength / 2, int(np.sqrt(NbRays)))
    RayList.append(mray.Ray(np.array([0, 0, 0]), np.array([0, 0, 1]), Number=0))
    num = 1
    for i in x:
        for j in y:
            if abs(x) > 1e-4 and abs(y) > 1e-4:
                RayList.append(mray.Ray(np.array([i, j, 0]), np.array([0, 0, 1]), Number=num, Wavelength=Wavelength))
                num += 1
    RayList = mgeo.RotationRayList(RayList, np.array([0, 0, 1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, Centre)
    return RayList


# %%
# def GaussianIntensity(r,z,I0, Wavelength, Divergence):

#     w0 = Wavelength / (np.pi * Divergence)
#     zR = np.pi * w0**2 / Wavelength
#     wz = w0 * np.sqrt(1+(z/zR)**2)
#     return I0 * np.exp(-2 * (r/wz)**2) / (1+(z/zR)**2)


def ApplyGaussianIntensityToRayList(RayList, IntensityFraction=1 / np.e**2):
    """
    Apply a Gaussain intensity profile to a ray list, with intensity =1 at the center and falling
    to 'IntensityFraction' at the edge (arbitrary units).

    Parameters
    ----------
        RayList : list of Ray-objects
            Coordinate vector of the center of the source disk (in mm).

        IntensityFraction : float, optional
            Relative intensity in (0,1), reached at the edge of the ray bundle. Default is 1/e^2.

    """
    if IntensityFraction >= 1 or IntensityFraction <= 0:
        # raise ValueError('IntensityFraction should be between 0 and 1!')
        print(
            "When applying a Gaussian intensity profile to a ray list, the IntensityFraction should be between 0 and 1! I'm setting it to 1/e^2."
        )
        IntensityFraction = 1 / np.e**2

    Axis = mp.FindCentralRay(RayList).vector
    Divergence = 0
    for l in RayList:
        # find the largest angle with the central ray out of all the rays in the list
        Divergence = max(Divergence, mgeo.AngleBetweenTwoVectors(Axis, l.vector))

    if (
        Divergence > 1e-12
    ):  # we're dealing with a point source and not a plane wave, so we can apply Gaussian intensity profile as function of angle
        for k in RayList:
            Angle = mgeo.AngleBetweenTwoVectors(Axis, k.vector)
            Intensity = np.exp(-2 * (np.tan(Angle) / Divergence) ** 2 * -0.5 * np.log(IntensityFraction))
            k.intensity = Intensity
    else:  # otherwise we apply it as function of ray distance from the Axis
        MaxDist = 0
        for l in RayList:
            MaxDist = max(MaxDist, np.linalg.norm(l.point))
        for k in RayList:
            Intensity = np.exp(-2 * (np.linalg.norm(k.point) / MaxDist) ** 2 * -0.5 * np.log(IntensityFraction))
            k.intensity = Intensity

    return RayList
