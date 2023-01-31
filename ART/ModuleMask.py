# -*- coding: utf-8 -*-
"""
Created in 2021

@author: Stefan Haessler
"""
#%% Modules

import numpy as np
import ART.ModuleGeometry as mgeo


#%%############################################################################
class Mask:
    """ a mask that lets all rays that DON'T hit it pass and deletes the rest """
    def __init__(self, Support):
        self.type = 'Mask'
        self.support = Support
        
    
    def get_intersection(self, Ray):
        """  Return the intersection point between Ray and the xy-plane,
        where rays do NOT hit the mask's support """
        t = -Ray.point[2] / Ray.vector[2]
        I = Ray.vector*t + Ray.point
        if t > 0 and not self.support._IncludeSupport(I):
            PointIntersection = I
        else: 
            PointIntersection = None
        
        return PointIntersection
    
    def get_normal(self, Point):
        Normal = np.array([0,0,1])
        return Normal
    
    def get_centre(self):
        return np.array([0,0,0])
    
    def get_grid3D(self,NbPoint):
        ListCoordXYZ = []
        ListCoordXY = self.support._get_grid(NbPoint)
        for k in ListCoordXY:
            z = 0
            ListCoordXYZ.append(np.array([k[0],k[1],z]))
        return ListCoordXYZ
 
    
    
#%%############################################################################
def TransmitMaskRay(Mask, PointMask, Ray):
    """ Returns the transmitted ray """   
    PointRay = Ray.point
    VectorRay = Ray.vector
    NormalMask = Mask.get_normal(PointMask)
        
    AngleIncidence = mgeo.AngleBetweenTwoVectors(VectorRay, NormalMask)
    Path = np.linalg.norm(PointMask - PointRay) + Ray.path
    
    RayTransmitted = Ray.copy_ray()
    RayTransmitted.point = PointMask
    RayTransmitted.vector = VectorRay
    RayTransmitted.incidence = AngleIncidence
    RayTransmitted.path = Path
    
    return RayTransmitted

#%%
def TransmitMaskRayList(Mask, RayList):
    
    ListRayTransmitted = []
    for Ray in RayList:
        PointMask = Mask.get_intersection(Ray)

        if PointMask is not None: 
            RayTransmitted = TransmitMaskRay(Mask, PointMask, Ray)
            ListRayTransmitted.append(RayTransmitted)

    return ListRayTransmitted