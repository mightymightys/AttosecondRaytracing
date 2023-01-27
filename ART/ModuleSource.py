#%% Modules

import numpy as np
import ART.ModuleOpticalRay as mray
import ART.ModuleGeometry as mgeo

#%%

def Cone(Angle, NbRays):
    """
    Return a list of rays filling a cone according to Vogel's spiral.
    Angle is the cone's angle. NbRays is the number of rays.
    """
    Height = 1
    Radius = Height * np.tan(Angle)
    
    MatrixXY = mgeo.SpiralVogel(NbRays, Radius)
    
    RayList = []
    
    for k in range(NbRays):
        x = MatrixXY[k,0]
        y = MatrixXY[k,1]
        RayList.append(mray.Ray(np.array([0,0,0]), np.array([x,y,Height]), Number=k))
  
    return RayList


def ExtendedSource(S, Axis, Diameter, Divergence, NbRays):
    """
    Return a list of rays simulating rays from an extended array of point sources, distributed over a disk with Diameter.
    S is the coordinate vector of the source, Axis is the axis vector of the cone (main direction of the rays), Divergence is the divergence of the source (similar to gaussian beam)
    """
    NbPointSources = max(5, int(250*Diameter))
    NbPointSources = min(NbPointSources, 100)
    
    MatrixXY = mgeo.SpiralVogel(NbPointSources, Diameter/2) #the positions of the point sources
    
    NbRaysPerPointSource = max(100, int(NbRays/NbPointSources))
    RayList = []
    PointSourceRayList = Cone(Divergence, NbRaysPerPointSource)
    for k in range(NbPointSources):
        ShiftedPointSourceRayList = mgeo.TranslationRayList(PointSourceRayList, [MatrixXY[k,0],MatrixXY[k,1],0])
        for l in range(NbRaysPerPointSource):
            ShiftedPointSourceRayList[l]._number += k*NbRaysPerPointSource #we allow ourselves exceptionally to modify the "internal" _number attribute
        RayList.extend(ShiftedPointSourceRayList)    
    
    RayList = mgeo.RotationRayList(RayList, np.array([0,0,1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, S)
    return RayList


def PointSource(S, Axis, Divergence, NbRays):
    """
    Return a list of rays simulating rays from a point source.
    S is the coordinate vector of the source, Axis is the axis vector of the cone (main direction of the rays), Divergence is the divergence of the source (similar to gaussian beam)
    """
    RayList = Cone(Divergence, NbRays)
    RayList = mgeo.RotationRayList(RayList, np.array([0,0,1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, S)
    return RayList

def GaussianIntensity(r,z,I0, Wavelength, Divergence):

    w0 = Wavelength / (np.pi * Divergence)
    zR = np.pi * w0**2 / Wavelength
    wz = w0 * np.sqrt(1+(z/zR)**2)
    return I0 * np.exp(-2 * (r/wz)**2) / (1+(z/zR)**2)


def ApplyGaussianIntensityToRayList(RayList, IntensityFraction):
    """
    Apply a gaussain intensity profile to a ray list (arbitrary units)
    """
    if IntensityFraction >=1 or IntensityFraction <=0:
        raise ValueError('IntensityFraction should be between 0 and 1!')
    Axis = RayList[0].vector
    Divergence = 0
    for l in RayList:
        Divergence = max(Divergence,mgeo.AngleBetweenTwoVectors(Axis, l.vector))
        
    if Divergence >1e-12: #we're dealing with a point source and not a plane wave, so we can apply Gaussian intensity  profile as function of angle
        for k in RayList:
            Angle = mgeo.AngleBetweenTwoVectors(Axis,k.vector)
            Intensity = np.exp(-2*(np.tan(Angle)/Divergence)**2 *-0.5*np.log(IntensityFraction)) 
            k.intensity = Intensity
    else:  #otherwise we apply it as function of ray distance from the Axis
        MaxDist = 0
        for l in RayList:
            MaxDist = max(MaxDist,np.linalg.norm(l.point))
        for k in RayList: 
            Intensity = np.exp(-2*(np.linalg.norm(k.point)/MaxDist)**2 *-0.5*np.log(IntensityFraction)) 
            k.intensity = Intensity
    return RayList
    
    
#%%
def PlaneWaveDisk(Centre, Axis, Radius, NbRays):
    MatrixXY = mgeo.SpiralVogel(NbRays, Radius)
    RayList = []
    num=0
    for k in range(NbRays-1):
        x = MatrixXY[k,0]
        y = MatrixXY[k,1]
        RayList.append(mray.Ray(np.array([x,y,0]),np.array([0,0,1]), Number=num))
        num = num + 1
    RayList = mgeo.RotationRayList(RayList, np.array([0,0,1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, Centre)
    return RayList

#%%
def PlaneWaveSquare(Centre, Axis, SideLength, NbRays):
    RayList = []
    x = np.linspace(-SideLength/2,SideLength/2,int(np.sqrt(NbRays)))
    y = np.linspace(-SideLength/2,SideLength/2,int(np.sqrt(NbRays)))
    RayList.append(mray.Ray(np.array([0,0,0]),np.array([0,0,1]), Number=0))
    num = 1
    for i in x:
        for j in y:
            if abs(x) >1e-4 and abs(y)>1e-4:
                RayList.append(mray.Ray(np.array([i,j,0]),np.array([0,0,1]), Number=num))
                num = num + 1
    RayList = mgeo.RotationRayList(RayList, np.array([0,0,1]), Axis)
    RayList = mgeo.TranslationRayList(RayList, Centre)
    return RayList

