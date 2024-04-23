"""
Provides a class which represents a virtual detector —a plane in 3D space— with methods
to analyze the impact points of a ray bundle intercepting it.

This module is used by 'ARTmain' to set up a detector according to the values in the *DetectorOptions*-dictionary filled in the CONFIG-scripts..

![Illustration of the detector plane.](detector.svg)



Created in Sept 2020

@author: Stefan Haessler
"""
# %% Modules

import numpy as np
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo

LightSpeed = 299792458000


# %%
class Detector:
    """
    A detector is simply a plane in 3D-lab-space.
    It is defined by a point Detector.centre and its
    normal vector Detector.normal (thought to point towards the incoming rays).
    These are optional arguments when creating a Detector-instance, because they can
    later be set automatically by the method Detector.autoplace().
    The Detector has a reference point "RefPoint", with respect to which a distance
    can be returned via the method Detector.get distance.

    Attributes
    ----------
        centre : np.ndarray
            3D Point in the Detector plane.

        normal : np.ndarray
           3D Normal vector on the Detector plane.

        refpoint : np.ndarray
            3D reference point from which to measure the Detector distance.
    """

    def __init__(self, RefPoint: np.ndarray, Centre=None, Normal=None):
        """
        Parameters
        ----------
            RefPoint : np.ndarray
                3D Reference point from which to measure the Detector distance.

            Centre : np.ndarray, optional
               3D Point in the Detector plane.

            Normal : np.ndarray, optional
               3D Normal vector on the Detector plane.
        """
        self.centre = Centre
        self.normal = Normal
        self.refpoint = RefPoint

    @property
    def centre(self):
        return self._centre

    # a setter function
    @centre.setter
    def centre(self, Centre):
        if type(Centre) == np.ndarray and Centre.shape == (3,) or Centre is None:
            self._centre = Centre
        else:
            raise TypeError("Detector Centre must be a 3D-vector, given as numpy.ndarray of shape (3,).")

    @property
    def normal(self):
        return self._normal

    @normal.setter
    def normal(self, Normal):
        if type(Normal) == np.ndarray and Normal.shape == (3,) and np.linalg.norm(Normal) > 0:
            self._normal = Normal / np.linalg.norm(Normal)
        elif Normal is None:
            self._normal = Normal
        else:
            raise TypeError("Detector Normal must be a 3D-vector of norm >0, given as numpy.ndarray of shape (3,).")

    @property
    def refpoint(self):
        return self._refpoint

    # a setter function
    @refpoint.setter
    def refpoint(self, RefPoint):
        if type(RefPoint) == np.ndarray and RefPoint.shape == (3,):
            self._refpoint = RefPoint
        else:
            raise TypeError("Detector RefPoint must a 3D-vector, given as numpy.ndarray of shape (3,).")

    # %% METHODS FOR PLACING THE DETECTOR ##################################################################

    def copy_detector(self):
        """
        Returns a new Detector object with the same properties.
        """
        return Detector(self.refpoint, self.centre, self.normal)

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
        CentralRay = mp.FindCentralRay(RayList)
        if CentralRay is None:
            DetectorNormal = np.array([0, 0, 0])
            CentralPoint = np.array([0, 0, 0])
            for k in RayList:
                DetectorNormal = DetectorNormal + k.vector
                CentralPoint = CentralPoint + k.point
            DetectorNormal = -DetectorNormal / len(RayList)
            CentralPoint = CentralPoint / len(RayList)
        else:
            DetectorNormal = -CentralRay.vector
            CentralPoint = CentralRay.point

        self.normal = DetectorNormal
        self.centre = CentralPoint - DetectorNormal * DistanceDetector
        self.refpoint = CentralPoint

    def get_distance(self):
        """
        Return distance of the Detector-plane from its reference point Detector.refpoint.
        """
        return np.linalg.norm(
            self.refpoint - mgeo.IntersectionLinePlane(self.refpoint, -self.normal, self.centre, self.normal)
        )

    def shiftToDistance(self, NewDistance: float):
        """
        Shift the Detector to the distance NewDistance from its reference point Detector.refpoint.

        Parameters
        ----------
            NewDistance : float
                The distance (absolute) to which to shift the Detector.
        """
        if type(NewDistance) == int or type(NewDistance) == float or type(NewDistance) == np.float64:
            shiftby = NewDistance - self.get_distance()  # by that much in the self.normal-direction
        else:
            raise TypeError("The new Detector Distance must be int or float.")

        self._centre += -shiftby * self.normal

    def shiftByDistance(self, Shift: float):
        """
        Shift the Detector *by* the distance Shift along its normal vector Detector.normal.

        Parameters
        ----------
            Shift : float
                The distance (relative) *by* which to shift the Detector.
        """
        if type(Shift) == int or type(Shift) == float or type(Shift) == np.float64:
            self._centre = self.centre - Shift * self.normal  # directly set _centre, i.e. circumvent setter
            # because we already checked shift here, and centre
            # and normal have been checked before
        else:
            raise TypeError("The Detector Distance Shift must be int or float.")

    def _iscomplete(self):
        """
        Checks whether the Detector plane is defined by a centre-point and normal-vector.
        """
        if self.centre is None or self.normal is None:
            raise TypeError("The detector has no centre and normal vectors defined yet.")
            return False
        else:
            return True

    # %% METHODS TO GET THE DETECTOR REPONSE ##################################################################

    def get_PointList3D(self, RayList) -> list[np.ndarray]:
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
        if self._iscomplete():
            ListPointDetector3D = []
            append = ListPointDetector3D.append
            for k in RayList:
                append(mgeo.IntersectionLinePlane(k.point, k.vector, self.centre, self.normal))
            return ListPointDetector3D

    def get_PointList2D(self, RayList) -> list[np.ndarray]:
        """
        Returns the list of 2D-points in the detector plane, with the origin at Detector.centre.

        Parameters
        ----------
            RayList : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            ListPointDetector2D : list[np.ndarray of shape (2,)]
        """
        if self._iscomplete():
            ListPointDetector3D = self.get_PointList3D(RayList)

            ListPointDetector3D = [
                k - self.centre for k in ListPointDetector3D
            ]  # Equivalent to taking the detector's central point as the origin
            ListPointDetector3D = mgeo.RotationPointList(ListPointDetector3D, self.normal, np.array([0, 0, 1]))

            ListPointDetector2D = [k[0:2] for k in ListPointDetector3D]
            return ListPointDetector2D

    def get_PointList2DCentre(self, RayList) -> list[np.ndarray]:
        """
        Returns the list of 2D-points in the detector plane, with the origin centered in the point cloud.

        Parameters
        ----------
            RayLis* : list[Ray]
                A list of objects of the ModuleOpticalRay.Ray-class.

        Returns
        ----------
            ListPointDetector2DCentre : list[np.ndarray of shape (2,)]
        """
        if self._iscomplete():
            ListPointDetector2D = self.get_PointList2D(RayList)
            ListPointDetector2DCentre = mgeo.CentrePointList(ListPointDetector2D)
            return ListPointDetector2DCentre

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
        if self._iscomplete():
            DetectorPointList3D = self.get_PointList3D(RayList)
            OpticalPathList = []
            for k in range(len(RayList)):
                OpticalPathList.append(np.linalg.norm(RayList[k].point - DetectorPointList3D[k]) + np.sum(RayList[k].path))

            MeanPath = np.mean(OpticalPathList)
            DelayList = [(k - MeanPath) / LightSpeed * 1e15 for k in OpticalPathList]  # in fs, c in mm/s
            return DelayList
