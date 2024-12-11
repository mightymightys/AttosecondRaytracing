"""
Provides a class which represents a virtual detector —a plane in 3D space— with methods
to analyze the impact points of a ray bundle intercepting it.

![Illustration of the detector plane.](detector.svg)



Created in Sept 2020

@author: Stefan Haessler + André Kalouguine
"""
# %% Modules

import numpy as np
import quaternion
from abc import ABC, abstractmethod

import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleOpticalRay as mray
import ARTcore.ModuleGeometry as mgeo
from ARTcore.ModuleGeometry import Point, Vector, Origin

# This is to allow multiple constructors based on the type of the first argument
from functools import singledispatchmethod
import logging
from copy import copy
import time


logger = logging.getLogger(__name__)

LightSpeed = 299792458000

# %% Abstract class
class Detector(ABC):
    """
    Abstract class describing a detector and the properties/methods it has to have. 
    It implements the required components for the detector to behave as a 3D object.
    Non-abstract subclasses of `Detector` should implement the remaining
    functiuons such as the detection of rays
    """
    @property
    def r(self):
        """
        Return the offset of the optical element from its reference frame to the lab reference frame.
        Not that this is **NOT** the position of the optical element's center point. Rather it is the
        offset of the reference frame of the optical element from the lab reference frame.
        """
        return self._r

    @r.setter
    def r(self, NewPosition):
        """
        Set the offset of the optical element from its reference frame to the lab reference frame.
        """
        if isinstance(NewPosition, mgeo.Point) and len(NewPosition) == 3:
            self._r = NewPosition
        else:
            raise TypeError("Position must be a 3D mgeo.Point.")
    @property
    def q(self):
        """
        Return the orientation of the optical element.
        The orientation is stored as a unit quaternion representing the rotation from the optic's coordinate frame to the lab frame.
        """
        return self._q
    
    @q.setter
    def q(self, NewOrientation):
        """
        Set the orientation of the optical element.
        This function normalizes the input quaternion before storing it.
        If the input is not a quaternion, raise a TypeError.
        """
        if isinstance(NewOrientation, np.quaternion):
            self._q = NewOrientation.normalized()
        else:
            raise TypeError("Orientation must be a quaternion.")

    @property
    def position(self):
        """
        Return the position of the basepoint of the optical element. Often it is the same as the optical centre.
        This position is the one around which all rotations are performed.
        """
        return self.r0 + self.r
    
    @property
    def orientation(self):
        """
        Return the orientation of the optical element.
        The utility of this method is unclear.
        """
        return self.q
    
    @property
    def basis(self):
        return self.r0, self.r, self.q
    
    def __init__(self):
        self._r = mgeo.Vector([0.0, 0.0, 0.0])
        self._q = np.quaternion(1, 0, 0, 0)
        self.r0 = mgeo.Point([0.0, 0.0, 0.0])

    @abstractmethod
    def centre(self):
        """
        (0,0) point of the detector.
        """
        pass

    @abstractmethod
    def refpoint(self):
        """
        Reference point from which to measure the distance of the detector. 
        Usually the center of the previous optical element.
        """
        pass

    @property
    def distance(self):
        """
        Return distance of the Detector from its reference point Detector.refpoint.
        """
        return (self.refpoint - self.centre).norm

    @distance.setter
    def distance(self, NewDistance: float):
        """
        Shift the Detector to the distance NewDistance from its reference point Detector.refpoint.

        Parameters
        ----------
            NewDistance : number
                The distance (absolute) to which to shift the Detector.
        """
        vector = (self.centre - self.refpoint).normalized
        self.centre = self.refpoint + vector * NewDistance

    @abstractmethod
    def get_2D_points(self, RayList):
        pass

    @abstractmethod
    def get_3D_points(self, RayList):
        pass

    @abstractmethod
    def __copy__(self):
        pass


# %% Infinite plane detector class
class InfiniteDetector(Detector):
    """
    Simple infinite plane.
    Beyond being just a solid object, the 

    Attributes
    ----------
        centre : np.ndarray
            3D Point in the Detector plane.

        refpoint : np.ndarray
            3D reference point from which to measure the Detector distance.
    """
    def __init__(self):
        super().__init__()
        self._centre = mgeo.Origin
        self._refpoint = mgeo.Origin
    
    def __copy__(self):
        """
        Returns a new Detector object with the same properties.
        """
        result = InfiniteDetector()
        result._centre = copy(self._centre)
        result._refpoint = copy(self._refpoint)
        result._r = copy(self._r)
        result._q = copy(self._q)
        result.r0 = copy(self.r0)
        return result
    
    @property
    def normal(self):
        return mgeo.Vector([0,0,1]).from_basis(*self.basis)

    @property
    def centre(self):
        return self._centre

    # a setter function
    @centre.setter
    def centre(self, Centre):
        if isinstance(Centre, np.ndarray) and Centre.shape == (3,):
            self._centre = mgeo.Point(Centre)
            self.r = self._centre # Simply because r0 is 0,0,0 anyways
        else:
            raise TypeError("Detector Centre must be a 3D-vector, given as numpy.ndarray of shape (3,).")
    
    @property
    def refpoint(self):
        return self._refpoint

    # a setter function
    @refpoint.setter
    def refpoint(self, RefPoint):
        if isinstance(RefPoint, np.ndarray) and RefPoint.shape == (3,):
            self._refpoint = mgeo.Point(RefPoint)
        else:
            raise TypeError("Detector RefPoint must a 3D-Point")

    # %% Detector placement methods

    @property
    def distance(self):
        """
        Return distance of the Detector from its reference point Detector.refpoint.
        """
        return (self.refpoint - self.centre).norm
    
    @distance.setter
    def distance(self, NewDistance: float):
        """
        Shift the Detector to the distance NewDistance from its reference point Detector.refpoint.

        Parameters
        ----------
            NewDistance : number
                The distance (absolute) to which to shift the Detector.
        """
        vector = (self.centre - self.refpoint).normalized()
        self.centre = self.refpoint + vector * NewDistance
    
    def set_distance(self,x):
        self.distance = x

    def autoplace(self, RayList, DistanceDetector: float):
        """
        Automatically place and orient the detector such that it is normal to the central ray
        of the supplied RayList, at the distance DistanceDetector the origin point of that central ray.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

            DistanceDetector : float
                The distance at which to place the Detector.
        """
        #CentralRay = mp.FindCentralRay(RayList)
        CentralRay = RayList[0]
        if CentralRay is None:
            logger.warning(f"Could not find central ray! The list of rays has a length of {len(RayList)}")
            CentralPoint = mgeo.Origin
            for k in RayList:
                CentralPoint = CentralPoint + (k.point - mgeo.Origin)
            CentralPoint = CentralPoint / len(RayList)
        else:
            logger.debug(f"Found central ray, using it to position detector: \n{CentralRay}")
            CentralPoint = CentralRay.point

        self.centre = CentralPoint + CentralRay.vector * DistanceDetector
        self.refpoint = CentralPoint
        self.q = mgeo.QRotationVector2Vector(mgeo.Vector([0,0,1]), -CentralRay.vector)

    # %% Detector optimisation methods
    def test_callback_distances(self, RayList, distances, callback,
                                provide_points=False,
                                detector_reference_frame=False):
        LocalRayList = [k.to_basis(*self.basis) for k in RayList]
        N = len(distances)
        # Calculate the impact points
        Points = mgeo.IntersectionRayListZPlane(LocalRayList, distances)
        values = []
        for k in range(N):
            Points[k] = Points[k]._add_dimension() + mgeo.Point([0,0,distances[k]])
            values += [callback(RayList, Points[k], self.basis)]
        return values

    def optimise_distance(self, RayList, Range,  callback, 
                          detector_reference_frame=False,
                          provide_points=False,
                          maxiter=5, 
                          tol=1e-6,
                          splitting=50,
                          Nrays=1000,
                          callback_iteration=None
        ):
        """
        Optimise the position of the detector within the provided range, trying to 
        maximise some value calculated by the callback function.

        The callback function receives a list of rays that already hit the detector.
        If requested, the rays can be provided in the reference frame of the detector.

        The callback function should return a single value that is to be minimised.

        The function splits the range into `splitting` points and uses the `IntersectionRayListZPlane`
        function to calculate the impact points for all of these.
        If `provide_points` is set to True, the function will pass on the result to the callback function.
        Otherwise it will generate rays for the callback function, including calculating the paths.

        It will then find the best value and redo the iteration around that value.
        Repeat until either the maximum number of iterations is reached or the tolerance is met.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

            Range : tuple
                The range of distances to consider for the detector placement.

            callback : function
                The function to be minimised. It should take a list of rays and return a single value.

            detector_reference_frame : bool
                If True, the rays are provided in the reference frame of the detector.

            provide_points : bool
                If True, the callback function will receive the impact points directly.

            maxiter : int
                The maximum number of iterations to perform.

            tol : float
                The tolerance to reach before stopping the iteration.

            splitting : int
                The number of points to split the range into.
        
        Returns
        ----------
            The optimal distance of the detector.
        """
        Range = [i-self.distance for i in Range]
        if Nrays<len(RayList):
            RayList = np.random.choice(RayList, Nrays, replace=False)
        previous_best = None
        values = [0]*splitting
        for i in range(maxiter):
            if callback_iteration is not None:
                callback_iteration(i, Range)
            # Split the range into `splitting` points
            distances = np.linspace(*Range, splitting)
            values = self.test_callback_distances(RayList, distances, callback, 
                                                  provide_points=provide_points, 
                                                  detector_reference_frame=detector_reference_frame)
            best = np.argmin(values)
            # Update the range
            if best == 0:
                Range = (distances[0], distances[2])
            elif best == splitting - 1:
                Range = (distances[-3], distances[-1])
            else:
                Range = (distances[best - 1], distances[best + 1])
            # Check if the tolerance is met
            if i>0:
                if np.abs(values[best] - previous_best) < tol:
                    break
            previous_best = values[best]
        self.distance -= distances[best]

    def _spot_size(self, RayList, Points, basis):
        """
        Returns focal spot size.

        Parameters
        ----------
            Points : mgeo.PointArray
                The points where the rays hit the detector.
        Returns
        ----------
            float
        """
        center = mgeo.Point(np.mean(Points[:,:2], axis=0))
        distances = (Points[:,:2] - center).norm
        return np.std(distances)
    
    def _delay_std(self, RayList, Points, basis):
        """
        Returns the standard deviation of the delays of the rays

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            float
        """
        #paths = np.array([np.sum(k.path) for k in RayList])
        paths = np.sum(np.array([k.path for k in RayList], dtype=float), axis=1)
        StartingPoints = mgeo.PointArray([k.point for k in RayList])
        XYZ = Points.from_basis(*basis)
        LastDistances = (XYZ - StartingPoints).norm
        return float(np.std(paths+LastDistances))

    def _over_intensity(self, RayList, Points, Basis):
        """
        Calculates spot_size * delay_std for the given RayList.
        """
        spot_size = self._spot_size(RayList, Points, Basis)
        delay_std = self._delay_std(RayList, Points, Basis)
        return spot_size * delay_std

    # %% Detector response methods 
    def get_3D_points(self, RayList) -> list[np.ndarray]:
        """
        Returns the list of 3D-points in lab-space where the rays in the
        list RayList hit the detector plane.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            ListPointDetector3D : list[np.ndarray of shape (3,)]
        """
        return self.get_2D_points(RayList)[0]._add_dimension().from_basis(*self.basis)

    def get_2D_points(self, RayList) -> list[np.ndarray]:
        """
        Returns the list of 2D-points in the detector plane, with the origin at Detector.centre.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class of length N

        Returns
        ----------
            XY : np.ndarray of shape (N,2)
        """
        return mgeo.IntersectionRayListZPlane([r.to_basis(*self.basis) for r in RayList])

    def get_centre_2D_points(self, RayList) -> list[np.ndarray]:
        """
        Returns the center of the 2D array of points.

        Parameters
        ----------
            RayList* : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            ListPointDetector2DCentre : list[np.ndarray of shape (2,)]
        """
        return np.mean(self.get_2D_points(RayList), axis=0)

    def get_Delays(self, RayList) -> list[float]:
        """
        Returns the list of delays of the rays of RayList as they hit the detector plane,
        calculated as their total path length divided by the speed of light.
        The delays are relative to the mean “travel time”, i.e. the mean path length of RayList
        divided by the speed of light.
        The list index corresponds to that of the RayList.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            DelayList : list[float]
        """
        XYZ = self.get_2D_points(RayList)[0]._add_dimension().from_basis(*self.basis)
        StartingPoints = mgeo.PointArray([k.point for k in RayList])
        LastDistances = (XYZ - StartingPoints).norm
        PreviousDistances = np.array([np.array(k.path) for k in RayList])
        Distances = np.hstack((PreviousDistances, LastDistances[:, np.newaxis]))
        TotalPaths = np.sum(Distances, axis=1)
        MeanPath = np.mean(TotalPaths)
        return list((TotalPaths-MeanPath) / LightSpeed * 1e15) # in fs