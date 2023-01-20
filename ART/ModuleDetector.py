"""
Created in Sept 2020

@author: Stefan Haessler
"""
#%% Modules

import numpy as np
import ART.ModuleProcessing as mp
import ART.ModuleGeometry as mgeo

LightSpeed = 299792458000

#%%
class Detector:
    """ 
    A detector is simply a plane in 3D placed in the optical chain.
    It is defined by a point "DetectorCentre" and its normal vector "DetectorNormal".
    The normal is thought to point towards the incoming rays.
    It has a  reference point "RefPoint", with respect to which a distance can be returned.
    This can e.g. be the centre of the preceding optical element.
    """
 
    def __init__(self, RefPoint, Centre=None, Normal=None):
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
        else: raise TypeError('Detector Centre must be a 3D-vector, given as numpy.ndarray of shape (3,).')
        
         
    @property
    def normal(self): 
        return self._normal
       
    @normal.setter 
    def normal(self, Normal): 
        if type(Normal) == np.ndarray and Normal.shape == (3,) and np.linalg.norm(Normal) > 0 : 
              self._normal = Normal/np.linalg.norm(Normal)
        elif Normal is None:
              self._normal = Normal
        else: raise TypeError('Detector Normal must be a 3D-vector of norm >0, given as numpy.ndarray of shape (3,).')
        
        
    @property
    def refpoint(self): 
        return self._refpoint
       
    # a setter function 
    @refpoint.setter 
    def refpoint(self, RefPoint): 
        if type(RefPoint) == np.ndarray and RefPoint.shape == (3,) : 
              self._refpoint = RefPoint
        else: raise TypeError('Detector RefPoint must a 3D-vector, given as numpy.ndarray of shape (3,).')

        
#%% METHODS FOR PLACING THE DETECTOR ##################################################################
   
    def copy_detector(self):
        return Detector(self.refpoint, self.centre, self.normal)
     
    
    def autoplace(self, RayList, DistanceDetector):
        CentralRay = mp.FindCentralRay(RayList)
        if CentralRay is None:
            DetectorNormal = np.array([0,0,0])
            CentralPoint = np.array([0,0,0])
            for k in RayList:
                DetectorNormal = DetectorNormal + k.vector
                CentralPoint = CentralPoint + k.point
            DetectorNormal = -DetectorNormal/len(RayList)  
            CentralPoint = CentralPoint/len(RayList)  
        else:    
            DetectorNormal = -CentralRay.vector
            CentralPoint = CentralRay.point
        
        self.normal  = DetectorNormal
        self.centre  = CentralPoint - DetectorNormal*DistanceDetector


    def get_distance(self):
        return np.linalg.norm(self.refpoint-mgeo.IntersectionLinePlane(self.refpoint, -self.normal, self.centre, self.normal))       
    
    
    def shiftToDistance(self, NewDistance):
        if type(NewDistance) == int or type(NewDistance) == float or type(NewDistance) == np.float64: 
            shiftby = NewDistance - self.get_distance() #by that much in the self.normal-direction
        else: raise TypeError('The new Detector Distance must be int or float.')
        
        self._centre = self.centre - shiftby*self.normal

    def shiftByDistance(self, Shift):
        if type(Shift) == int or type(Shift) == float or type(Shift) == np.float64: 
            self._centre = self.centre - Shift*self.normal  #directly set _centre, i.e. circumvent setter
                                                            #because we already checked shift here, and centre
                                                            #and normal have been checked before
        else: raise TypeError('The Detector Distance Shift must be int or float.')
        

    def iscomplete(self):
        if self.centre is None or self.normal is None:
            raise TypeError('The detector has no centre and normal vectors defined yet.')
            return False
        else: return True
        
#%% METHODS TO GET THE DETECTOR REPONSE ##################################################################
        
    def get_PointList3D(self, RayList):
        if self.iscomplete():
            ListPointDetector3D = []
            for k in RayList:
                ListPointDetector3D.append(mgeo.IntersectionLinePlane(k.point, k.vector, self.centre, self.normal))
            return ListPointDetector3D
    
    def get_PointList2D(self, RayList):
        if self.iscomplete():
            ListPointDetector3D = self.get_PointList3D(RayList)
    
            ListPointDetector3D = [k - self.centre for k in ListPointDetector3D] # Equivalent to taking the detector's central point as the origin
            ListPointDetector3D = mgeo.RotationPointList(ListPointDetector3D, self.normal, np.array([0,0,1]))
    
            ListPointDetector2D = [k[0:2] for k in ListPointDetector3D]
            return ListPointDetector2D

    def get_PointList2DCentre(self, RayList):
        if self.iscomplete():
            ListPointDetector2D = self.get_PointList2D(RayList)
            ListPointDetector2DCentre = mp.CentrePointList(ListPointDetector2D)
            return ListPointDetector2DCentre


    def get_Delays(self, RayList):
        if self.iscomplete():
            DetectorPointList3D = self.get_PointList3D(RayList)
            OpticalPathList = []
            for k in range(len(RayList)):
                OpticalPathList.append(np.linalg.norm(RayList[k].point - DetectorPointList3D[k]) + RayList[k].path)
        
            MeanPath = np.mean(OpticalPathList)
            DelayList = [(k-MeanPath)/LightSpeed*1e15 for k in OpticalPathList] # in fs, c in mm/s
            return DelayList