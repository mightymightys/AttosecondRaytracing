"""
Contains processing functions used for the ray-tracing and automated (pre-)alignment of the optical chains.
Usually these don't need to be called by users of ART, but they may be useful.



Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
# %% Modules
import os
import pickle

# import gzip
import lzma
from time import perf_counter
from datetime import datetime
import copy
import numpy as np
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleMask as mmask
import ARTcore.ModuleSource as msource
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleOpticalChain as moc


# %%
def _singleOEPlacement(
    SourceProperties: dict,
    OpticsList: list,
    DistanceList: list[(int, float)],
    IncidenceAngleList: list[(int, float)],
    IncidencePlaneAngleList: list[(int, float)],
    Description: str,
):
    """
    Automatic placement and alignment of the optical elements for one optical chain.
    """
    Divergence = SourceProperties["Divergence"]
    SourceSize = SourceProperties["SourceSize"]
    RayNumber = SourceProperties["NumberRays"]
    Wavelength = SourceProperties["Wavelength"]

    IncidencePlaneAngleList = [
        np.deg2rad((i % 360)) for i in IncidencePlaneAngleList
    ]  # wrap angles into [0,360]deg and convert to radian
    IncidenceAngleList = [
        np.deg2rad((i % 360)) for i in IncidenceAngleList
    ]  # wrap angles into [0,360]deg and convert to radian

    # Launch source ray bundle from the origin [x,y,z]=[0,0,0] into the x-direction:
    SourcePosition = np.array([0, 0, 0])
    SourceDirection = np.array([1, 0, 0])
    if Divergence == 0:
        if SourceSize == 0:
            try:
                radius = 0.5 * min(OpticsList[0].support.dimX, OpticsList[0].support.dimY)  # for rect. support
            except AttributeError:
                radius = OpticsList[0].support.radius  # otherwise it must be a round support
        else:
            radius = SourceSize / 2
        SourceRayList = msource.PlaneWaveDisk(SourcePosition, SourceDirection, radius, RayNumber, Wavelength=Wavelength)
    else:
        if SourceSize == 0:
            SourceRayList = msource.PointSource(
                SourcePosition, SourceDirection, Divergence, RayNumber, Wavelength=Wavelength
            )
        else:
            SourceRayList = msource.ExtendedSource(
                SourcePosition, SourceDirection, SourceSize, Divergence, RayNumber, Wavelength=Wavelength
            )

    SourceRayList = msource.ApplyGaussianIntensityToRayList(
        SourceRayList, 1 / np.e**2
    )  # Intensity drops to 1/e^2 at edge of ray-bundle

    # Now successively create and align optical elements
    SourceRay = [
        mray.Ray(SourcePosition, SourceDirection)
    ]  # a single source ray in the central direction of the ray-bundle created above
    OpticalElements = []
    OpticalElementCentre = SourcePosition
    CentralVector = SourceDirection
    RotationAxis = np.array(
        [0, 1, 0]
    )  # Vector perpendicular to the incidence plane, i.e. initially the incidence plane is the x-z-plane.

    for k, Optic in enumerate(OpticsList):
        # for convex mirrors, rotated them by 180° them so we reflect from the "back side"
        if Optic.type == "SphericalCX Mirror" or Optic.type == "CylindricalCX Mirror":
            IncidenceAngleList[k] = np.pi - IncidenceAngleList[k]

        # shift OpticalElementCentre-point from that of the preceding OE, by the required distance along the central ray vector
        OpticalElementCentre = CentralVector * DistanceList[k] + OpticalElementCentre

        if abs(IncidencePlaneAngleList[k] - np.pi) < 1e-10:
            RotationAxis = -RotationAxis
        else:
            RotationAxis = mgeo.RotationAroundAxis(CentralVector, -IncidencePlaneAngleList[k], RotationAxis)

        OpticalElementNormal = mgeo.RotationAroundAxis(
            RotationAxis, -np.pi / 2 + IncidenceAngleList[k], np.cross(CentralVector, RotationAxis)
        )

        OpticalElementMajorAxis = np.cross(RotationAxis, OpticalElementNormal)

        Element = moe.OpticalElement(Optic, OpticalElementCentre, OpticalElementNormal, OpticalElementMajorAxis)

        OpticalElements.append(Element)
        auxChain = moc.OpticalChain(SourceRay, OpticalElements)

        if "Mirror" in Optic.type:
            OutRays = auxChain.get_output_rays()
            CentralVector = OutRays[-1][0].vector
        elif Optic.type == "Mask":
            # CentralVector remains unchanged
            # in our auxChain, but *not* in the OpticalElements-list, replace mask by completely tansparent "fake version"
            # to be sure that our CentralVector serving as alignment-guide always passes
            FakeMask = mmask.Mask(msupp.SupportRoundHole(Radius=100, RadiusHole=100, CenterHoleX=0, CenterHoleY=0))
            auxChain.optical_elements[-1] = moe.OpticalElement(
                FakeMask, OpticalElementCentre, OpticalElementNormal, OpticalElementMajorAxis
            )
        else:
            raise NameError("I don`t recognize the type of optical element " + OpticalElements[k].type.type + ".")

    return moc.OpticalChain(SourceRayList, OpticalElements, Description)


def OEPlacement(
    SourceProperties: dict,
    OpticsList: list,
    DistanceList: list[(int, float)],
    IncidenceAngleList: list[float],
    IncidencePlaneAngleList: list[float] = None,
    Description: str = ""
):
    """
    Automatically place optical elements in the "lab frame" according to given distances and incidence angles.
    Outputs an OpticalChain-object.

    The source is placed at the origin, and points into the x-direction. The optical elements then follow.

    As long as the angles in IncidencePlaneAngleList are 0, the incidence plane remains the x-z plane, i.e.
    the optical elements are rotated about the y-axis to set the desired incidence angle between OE-normal
    and the direction of incidence of the incoming ray-bundle. Otherwise, the incidence plane gets rotated.

    One of the elements of one of the lists 'DistanceList', 'IncidenceAngleList', or 'IncidencePlaneAngleList'
    can be a list or a numpy-array. In that case, a list of OpticalChain-objects is created which can be looped over.

    Parameters
    ----------
        SourceProperties : dict
            Dictionary with the keys "Divergence" :float, "SourceSize" :float,
            "NumberRays" :int, and "Wavelength" :float.
            This determines the source ray-bundle that is created.

        OpticsList : list[Optics-class-objects (Mirorrs or Masks)]
            List of Optics-objects.

        DistanceList : list[float]
            List of distances, in mm,  at which place the optics out of the
            OpticsList from the precending one. For the first optic, this is
            the distance from the light source.

            One of the elements can ne a list or a numpy-array instead of a float.
            In that case, a list of OpticalChain-objects is created which can be looped over.

        IncidenceAngleList : list[float]
            List of incidence angles, in degrees, measured between the central ray
            of the incident ray bundle and the normal on the optical element.

            One of the elements can ne a list or a numpy-array instead of a float.
            In that case, a list of OpticalChain-objects is created which can be looped over.

        IncidencePlaneAngleList : list[float], optional
            List of angles, in degrees, by which the incidence plane is rotated
            on the optical element with respect to that on the preceding one.
            Consequently, for the first element, this angle has no qualitative effect.

            One of the elements can ne a list or a numpy-array instead of a float.
            In that case, a list of OpticalChain-objects is created which can be looped over.

        Description : str, optional
            A string to describe the optical setup.

    Returns
    -------
        OpticalChainList : list[OpticalChain-object]

        or

        OpticalChain : OpticalChain-object

    """

    if IncidencePlaneAngleList is None:
        IncidencePlaneAngleList = np.zeros(len(OpticsList)).tolist()

    nest_indx_distance = _which_indeces(DistanceList)
    nest_indx_incidence = _which_indeces(IncidenceAngleList)
    nest_indx_incplane = _which_indeces(IncidencePlaneAngleList)

    total_nested = len(nest_indx_incidence + nest_indx_incplane + nest_indx_distance)

    if total_nested > 1:
        raise ValueError(
            "Only one element of one of the lists IncidenceAngleList, IncidencePlaneAngleList, or DistanceList can be a list or array itself. Otherwise things get too tangled..."
        )
    elif total_nested == 1:
        i = (nest_indx_incidence + nest_indx_incplane + nest_indx_distance)[0]
        loop_variable_name = OpticsList[i].type + "_idx_" + str(i)
        if nest_indx_incidence:  # if the list is not empty
            loop_variable_name += " incidence angle (deg)"
            loop_variable = copy.deepcopy(IncidenceAngleList[i])
            loop_list = IncidenceAngleList
        elif nest_indx_distance:  # if the list is not empty
            loop_variable_name += " distance (mm)"
            loop_variable = copy.deepcopy(DistanceList[i])
            loop_list = DistanceList
        elif nest_indx_incplane:  # if the list is not empty
            loop_variable_name += " incidence-plane angle rotation (deg)"
            loop_variable = copy.deepcopy(IncidencePlaneAngleList[i])
            loop_list = IncidencePlaneAngleList

        OpticalChainList = []
        for x in loop_variable:
            loop_list[i] = x
            ModifiedOpticalChain = _singleOEPlacement(
                SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description
            )
            ModifiedOpticalChain.loop_variable_name = loop_variable_name
            ModifiedOpticalChain.loop_variable_value = x
            OpticalChainList.append(ModifiedOpticalChain)

        return OpticalChainList

    elif total_nested == 0:
        OpticalChain = _singleOEPlacement(
            SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description
        )  # all simple
        
        return OpticalChain


# %%
def RayTracingCalculation(
    source_rays: list[mray.Ray], optical_elements: list[moe.OpticalElement], IgnoreDefects = True
) -> list[list[mray.Ray]]:
    """
    The actual ray-tracing calculation, starting from the list of 'source_rays',
    and propagating them from one optical element to the next in the order of the
    items of the list 'optical_elements'.


    Parameters
    ----------
        source_rays : list[Ray-objects]
            List of input rays.

        optical_elements : list[OpticalElement-objects]


    Returns
    -------
        output_rays : list[list[Ray-objects]]
            List of lists of rays, each item corresponding to the ray-bundle *after*
            the item with the same index in the 'optical_elements'-list.
    """
    ez = np.array([0, 0, 1])
    ex = np.array([1, 0, 0])
    output_rays = []

    for k in range(0, len(optical_elements)):
        if k == 0:
            RayList = source_rays
        else:
            RayList = output_rays[k - 1]

        # Mirror = optical_elements[k].type
        Position = optical_elements[k].position
        n = optical_elements[k].normal
        m = optical_elements[k].majoraxis

        # transform rays into mirror coordinate system
        RayList = mgeo.TranslationRayList(RayList, -Position)
        RayList = mgeo.RotationRayList(RayList, n, ez)
        # FIRST rotate m in the same way as you just did the ray list,
        # THEN rotate the rays again to match the NEW m with the x-axis !
        mPrime = mgeo.RotationPoint(m, n, ez)
        RayList = mgeo.RotationRayList(RayList, mPrime, ex)
        RayList = mgeo.TranslationRayList(RayList, optical_elements[k].type.get_centre())

        # optical element acts on the rays:
        if "Mirror" in optical_elements[k].type.type:
            RayList = mmirror.ReflectionMirrorRayList(optical_elements[k].type, RayList, IgnoreDefects = IgnoreDefects)
        elif optical_elements[k].type.type == "Mask":
            RayList = mmask.TransmitMaskRayList(optical_elements[k].type, RayList)
        else:
            raise NameError("I don`t recognize the type of optical element " + optical_elements[k].type.type + ".")

        # transform new rays back into "lab frame"
        RayList = mgeo.TranslationRayList(RayList, -optical_elements[k].type.get_centre())
        RayList = mgeo.RotationRayList(RayList, ex, mPrime)
        RayList = mgeo.RotationRayList(RayList, ez, n)
        RayList = mgeo.TranslationRayList(RayList, Position)

        output_rays.append(RayList)

    return output_rays


# %%
def _FindOptimalDistanceBIS(movingDetector, Amplitude, Step, RayList, OptFor, IntensityWeighted):
    ListSizeSpot = []
    ListDuration = []
    ListFitness = []
    if IntensityWeighted:
        Weights = [k.intensity for k in RayList]

    movingDetector.shiftByDistance(-Amplitude)
    n = int(2 * Amplitude / Step)
    for i in range(n):
        ListPointDetector2DCentre = movingDetector.get_PointList2DCentre(RayList)
        if OptFor in ["intensity", "spotsize"]:
            if IntensityWeighted:
                SpotSize = WeightedStandardDeviation(ListPointDetector2DCentre, Weights)
            else:
                SpotSize = StandardDeviation(ListPointDetector2DCentre)
            ListSizeSpot.append(SpotSize)

        if OptFor in ["intensity", "duration"]:
            DelayList = movingDetector.get_Delays(RayList)
            if IntensityWeighted:
                Duration = WeightedStandardDeviation(DelayList, Weights)
            else:
                Duration = StandardDeviation(DelayList)
            ListDuration.append(Duration)

        if OptFor == "intensity":
            Fitness = SpotSize**2 * Duration
        elif OptFor == "duration":
            Fitness = Duration
        elif OptFor == "spotsize":
            Fitness = SpotSize
        ListFitness.append(Fitness)

        movingDetector.shiftByDistance(Step)

    FitnessMin = min(ListFitness)
    ind = ListFitness.index(FitnessMin)
    if OptFor in ["intensity", "spotsize"]:
        OptSizeSpot = ListSizeSpot[ind]
    else:
        OptSizeSpot = np.nan
    if OptFor in ["intensity", "duration"]:
        OptDuration = ListDuration[ind]
    else:
        OptDuration = np.nan

    movingDetector.shiftByDistance(-(n - ind) * Step)

    return movingDetector, OptSizeSpot, OptDuration


def FindOptimalDistance(
    Detector,
    RayList: list[mray.Ray],
    OptFor="intensity",
    Amplitude: float = None,
    Precision: int = 3,
    IntensityWeighted=False,
    verbose=False,
):
    """
    Automatically finds the optimal the 'Detector'-distance for either
    maximum intensity, or minimal spotsize or duration.

    Parameters
    ----------
        Detector : Detector-object

        RayList : list[Ray-objects]

        OptFor : str, optional
            "spotsize": minimizes the standard-deviation spot size *d* on the detector.
            "duration": minimizes the standard-deviation *tau* of the ray-delays.
            "intensity": Maximizes 1/tau/d^2.
            Defaults to "intensity".

        Amplitude : float, optional
            The detector-distances within which the optimization is done will be
            the distance of 'Detector' +/- Amplitude in mm.

        Precision : int, optional
            Precision-parameter for the search algorithm. For Precision = n,
            it will iterate with search amplitudes decreasing until 10^{-n}.
            Defaults to 3.

        IntensityWeighted : bool, optional
            Whether to weigh the calculation of spotsize and/or duration by the ray-intensities.
            Defaults to False.

        verbose : bool, optional
            Whether to print results to the console.

    Returns
    -------
        OptDetector : Detector-object
            A new detector-instance with optimized distance.

        OptSpotSize : float
            Spotsize, calculated as standard deviation in mm of the spot-diagram
            on the optimized detector.

        OptDuration : float
            Duration, calculated as standard deviation in fs of the ray-delays
            on the optimized detector.
    """

    if OptFor not in ["intensity", "size", "duration"]:
        raise NameError(
            "I don`t recognize what you want to optimize the detector distance for. OptFor must be either 'intensity', 'size' or 'duration'."
        )

    FirstDistance = Detector.get_distance()
    ListPointDetector2DCentre = Detector.get_PointList2DCentre(RayList)
    SizeSpot = 2 * StandardDeviation(ListPointDetector2DCentre)
    NumericalAperture = ReturnNumericalAperture(RayList, 1)
    if Amplitude is None:
        Amplitude = min(4 * np.ceil(SizeSpot / np.tan(np.arcsin(NumericalAperture))), FirstDistance)
    Step = Amplitude/10
    #Step = Amplitude / 5  # pretty good in half the time ;-)

    if verbose:
        print(
            f"Searching optimal detector position for *{OptFor}* within [{FirstDistance-Amplitude:.3f}, {FirstDistance+Amplitude:.3f}] mm...",
            end="",
            flush=True,
        )
    movingDetector = Detector.copy_detector()
    for k in range(Precision + 1):
        movingDetector, OptSpotSize, OptDuration = _FindOptimalDistanceBIS(
            movingDetector, Amplitude * 0.1**k, Step * 0.1**k, RayList, OptFor, IntensityWeighted
        )

    if (
        not FirstDistance - Amplitude + 10**-Precision
        < movingDetector.get_distance()
        < FirstDistance + Amplitude - 10**-Precision
    ):
        print("There`s no minimum-size/duration focus in the searched range.")

    print(
        "\r\033[K", end="", flush=True
    )  # move to beginning of the line with \r and then delete the whole line with \033[K
    return movingDetector, OptSpotSize, OptDuration


# %%
def FindCentralRay(RayList: list[mray.Ray]):
    """
    Out of a the list of Ray-objects, RayList, determine the average direction
    vector and source point, and the return a Ray-object representing this 
    central ray of the RayList.
    
    Parameters
    ----------
        RayList : list of Ray-objects

    Returns
    -------
        CentralRay : Ray-object
    """
      
    CentralVector = np.mean( [x.vector for x in RayList], axis=0)
    CentralPoint = np.mean( [x.point for x in RayList], axis=0)
    
    return mray.Ray(CentralPoint, CentralVector)


def StandardDeviation(List: list[float, np.ndarray]) -> float:
    """
    Return the standard deviation of elements in the supplied list.
    If the elements are floats, it's simply the root-mean-square over the numbers.
    If the elements are vectors (represented as numpy-arrays), the root-mean-square
    of the distance of the points from their mean is calculated.

    Parameters
    ----------
        List : list of floats or of np.ndarray
            Numbers (in ART, typically delays)
            or 2D or 3D vectors (in ART typically space coordinates)

    Returns
    -------
        Std : float
    """
    if type(List[0]) in [int, float, np.float64]:
        return np.std(List)
    elif len(List[0]) > 1:
        return np.sqrt(np.var(List, axis=0).sum())
    else:
        raise ValueError("StandardDeviation expects a list of floats or numpy-arrays as input, but got something else.")


def WeightedStandardDeviation(List: list[float, np.ndarray], Weights: list[float]) -> float:
    """
    Return the weighted standard deviation of elements in the supplied list.
    If the elements are floats, it's simply the root-mean-square over the weighted numbers.
    If the elements are vectors (represented as numpy-arrays), the root-mean-square
    of the weighted distances of the points from their mean is calculated.

    Parameters
    ----------
        List : list of floats or of np.ndarray
            Numbers (in ART, typically delays)
            or 2D or 3D vectors (in ART typically space coordinates)

        Weights : list of floats, same length as List
            Parameters such as Ray.intensity to be used as weights

    Returns
    -------
        Std : float
    """
    average = np.average(List, axis=0, weights=Weights, returned=False)
    variance = np.average((List - average) ** 2, axis=0, weights=Weights, returned=False)
    return np.sqrt(variance.sum())


# %%
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
    CentralRay = FindCentralRay(RayList)
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


# %%
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


# %%
def _hash_list_of_objects(list):
    """returns a hash value for a list of objects, by just summing up all the individual hashes."""
    total_hash = 0
    for x in list:
        total_hash += hash(x)
    return total_hash


def _which_indeces(lst):
    """takes an input-list, and gives you a list of the indeces where the input-list contains another list or a numpy-array"""
    indexes = [i for i, x in enumerate(lst) if isinstance(x, (list, np.ndarray))]
    return indexes


# %%
def save_compressed(obj, filename: str = None):
    """Save (=pickle) an object 'obj' to a compressed file with name 'filename'."""
    if not type(filename) == str:
        filename = "kept_data_" + datetime.now().strftime("%Y-%m-%d-%Hh%M")

    i = 0
    while os.path.exists(filename + f"_{i}.xz"):
        i += 1
    filename = filename + f"_{i}"
    # with gzip.open(filename + '.gz', 'wb') as f:
    with lzma.open(filename + ".xz", "wb") as f:
        pickle.dump(obj, f)
    print("Saved results to " + filename + ".xz.")
    print("->To reload from disk do: kept_data = mp.load_compressed('" + filename + "')")


def load_compressed(filename: str):
    """Load (=unpickle) an object 'obj' from a compressed file with name 'filename'."""
    # with gzip.open(filename + '.gz', 'rb') as f:
    with lzma.open(filename + ".xz", "rb") as f:
        obj = pickle.load(f)
    return obj


# %%  tic-toc timer functions
_tstart_stack = []


def _tic():
    _tstart_stack.append(perf_counter())


def _toc(fmt="Elapsed: %s s"):
    print(fmt % (perf_counter() - _tstart_stack.pop()))
