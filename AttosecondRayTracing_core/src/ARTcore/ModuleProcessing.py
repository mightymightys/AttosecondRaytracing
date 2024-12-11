"""
Contains processing functions used for the ray-tracing and automated (pre-)alignment of the optical chains.
Usually these don't need to be called by users of ART, but they may be useful.



Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
# %% Modules
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleMask as mmask
import ARTcore.ModuleSource as msource
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleOpticalChain as moc

import os
import pickle

# import gzip
import lzma
from time import perf_counter
from datetime import datetime
import copy
import numpy as np
import quaternion
import logging

logger = logging.getLogger(__name__)


# %%
def singleOEPlacement(
    Optic: moe.OpticalElement,
    Distance: float,
    IncidenceAngle: float = 0,
    IncidencePlaneAngle: float = 0,
    InputRay: mray.Ray = None,
    AlignmentVector: str = "",
    PreviousIncidencePlane: mgeo.Vector = mgeo.Vector([0, 1, 0])
):
    """
    Automatic placement and alignment of a single optical element.
    The arguments are:
    - The optic to be placed.
    - The distance from the previous optic.
    - An incidence angle.
    - An incidence plane angle.
    - An input ray.
    - An alignment vector.

    The alignment procedure is as follows:
    - The optic is initially placed with its center at the source position (origin point of previous master ray). 
    - It's oriented such that the alignment vector is antiparallel to the master ray. By default, the alignment vector is the support_normal.
    - The majoraxis of the optic is aligned with the incidence plane. 
    - The optic is rotated around the incidence plane by the incidence angle.
    - The optic is rotated around the master ray by the incidence plane angle.
    - The optic is translated along the master ray by the distance from the previous optic.
    - The master ray is propagated through the optic without any blocking or diffraction effects and the output ray is used as the master ray for the next optic.
    """
    # Convert angles to radian and wrap to 2pi
    IncidencePlaneAngle = np.deg2rad(IncidencePlaneAngle) % (2 * np.pi)
    IncidenceAngle = np.deg2rad(IncidenceAngle) % (2 * np.pi)
    if InputRay is None:
        InputRay = mray.Ray(
            point=mgeo.Point([0, 0, 0]), 
            vector=mgeo.Vector([1, 0, 0]),
            path=(0.0,),
            number=0,
            wavelength=800e-9,
            incidence=0.0,
            intensity=1.0
    )
    OldOpticalElementCentre = InputRay.point
    MasterRayDirection = InputRay.vector.normalized()

    OpticalElementCentre = OldOpticalElementCentre + MasterRayDirection * Distance

    logger.debug(f"Old Optical Element Centre: {OldOpticalElementCentre}")
    logger.debug(f"Master Ray: {InputRay}")
    logger.debug(f"Optical Element Centre: {OpticalElementCentre}")
    # for convex mirrors, rotated them by 180° while keeping same incidence plane so we reflect from the "back side"
    if Optic.curvature == mmirror.Curvature.CONVEX:
        IncidenceAngle = np.pi - IncidenceAngle

    if hasattr(Optic, AlignmentVector):
        OpticVector = getattr(Optic, AlignmentVector)
    else:
        logger.warning(f"Optical Element {Optic} does not have an attribute {AlignmentVector}. Using support_normal instead.")
        OpticVector = Optic.support_normal_ref
    
    MajorAxis = Optic.majoraxis_ref

    IncidencePlane = PreviousIncidencePlane.rotate(mgeo.QRotationAroundAxis(InputRay.vector, -IncidencePlaneAngle))
    # We calculate a quaternion that will rotate OpticVector against the master ray and MajorAxis into the incidence plane
    # The convention is that the MajorAxis vector points right if seen from the source with IncidencePlane pointing up
    # To do that, we first calculate a vector MinorAxis in the reference frame of the optic that is orthogonal to both OpticVector and MajorAxis
    # Then we calculate a rotation matrix that will:
    #   - rotate MinorAxis into the IncidencePlane direction
    #   - rotate OpticVector against the master ray direction
    #   - rotate MajorAxis into the direction of the cross product of the two previous vectors
    MinorAxis = -np.cross(OpticVector, MajorAxis)
    qIncidenceAngle = mgeo.QRotationAroundAxis(IncidencePlane, -IncidenceAngle)
    q = mgeo.QRotationVectorPair2VectorPair(MinorAxis, IncidencePlane, OpticVector, -MasterRayDirection)
    q = qIncidenceAngle*q
    Optic.q = q
    logger.debug(f"Optic: {Optic}")
    logger.debug(f"OpticVector: {OpticVector}")
    logger.debug(f"Rotated OpticVector: {OpticVector.rotate(q)}")
    logger.debug(f"MasterRayDirection: {MasterRayDirection}")

    Optic.r = Origin + (OpticalElementCentre - Optic.centre - Optic.r)

    NextRay = Optic.propagate_raylist(mray.RayList.from_list([InputRay]), alignment=True)[0]
    
    return NextRay, IncidencePlane


def OEPlacement(
    OpticsList: list,
    InitialRay:mray.Ray = None,
    Source: msource.Source = None,
    InputIncidencePlane: mgeo.Vector = None
):
    """
    Automatic placement and alignment of the optical elements for one optical chain.
    Returns the output ray and the incidence plane of the last optical element.
    """
    if InitialRay is None:
        if Source is None:
            InputRay = InputRay = mray.Ray(
            point=mgeo.Point([0, 0, 0]), 
            vector=mgeo.Vector([1, 0, 0]),
            path=(0.0,),
            number=0,
            wavelength=800e-9,
            incidence=0.0,
            intensity=1.0
    )
        else:
            InputRay = Source.get_master_ray()
    else:
        InputRay = InitialRay.copy_ray()
    if InputIncidencePlane is None:
        InputIncidencePlane = mgeo.Vector([0, 1, 0])
    logger.debug(f"Initial Ray: {InputRay}")
    logger.debug(f"Initial Incidence Plane: {InputIncidencePlane}")
    assert np.linalg.norm(np.dot(InputRay.vector, InputIncidencePlane))<1e-6, "InputIncidencePlane is not orthogonal to InputRay.vector"
    PreviousIncidencePlane = InputIncidencePlane.copy()
    for i in range(len(OpticsList)):
        InputRay, PreviousIncidencePlane = singleOEPlacement(
            OpticsList[i]["OpticalElement"],
            OpticsList[i]["Distance"],
            OpticsList[i]["IncidenceAngle"],
            OpticsList[i]["IncidencePlaneAngle"],
            InputRay,
            OpticsList[i]["Alignment"] if "Alignment" in OpticsList[i] else "support_normal",
            PreviousIncidencePlane
        )
    return [i["OpticalElement"] for i in OpticsList]
    


# %%
def RayTracingCalculation(
    source_rays: mray.RayList, optical_elements: list[moe.OpticalElement], **kwargs
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
    output_rays = []

    for k in range(0, len(optical_elements)):
        if k == 0:
            RayList = source_rays
        else:
            RayList = output_rays[k - 1]
        RayList = optical_elements[k].propagate_raylist(RayList, **kwargs)
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
    
    return mray.Ray(mgeo.Point(CentralPoint), mgeo.Vector(CentralVector))


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
