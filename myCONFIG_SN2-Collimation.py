"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules

import numpy as np
#import addcopyfighandler  #if you want to be able to copy figure contents to the clipboard by ctrl-c, Only for windows though.

import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp


#%% User cell
verbose = True

FindBestSpot = False
FindBestDuration = False
FindOptDistance  = False

RayGraph = True
maxRaysToRender = 100

SpotDiagram = False
DelaySpotDiagram = True
IntensitySpotDiagram = False
IncidenceSpotDiagram = False

DelayGraph = False
IntensityGraph = False
IncidenceGraph = False
DrawAiryAndFourier = True

DelayMirrorProjection = True
IntensityMirrorProjection = False
IncidenceMirrorProjection = False


#%%
"""a telescope made from 2 spherical mirrors, and then a parabola"""

#%% optical chain
Divergence = 1.4e-3  
SourceSize = 340e-3   # in mm
Wavelength = 800e-6 
DeltaFT = 2
NumberRays = 2000

DistanceList = [3500, 1900, 550, 500, 1000]
AngleList = [-1, 46.5, 4, 3, 0.0]
IncidencePlaneAngleList = [0,90,-90,0,0]

Mirror1 = mmirror.MirrorSpherical(7000, msupp.SupportRound(12.5)) 
Mirror2 = mmirror.MirrorPlane(msupp.SupportRound(25)) 
Mirror3 = mmirror.MirrorSpherical(-1500, msupp.SupportRound(25)) 
Mirror4 = mmirror.MirrorSpherical(2500, msupp.SupportRound(25)) 

offAxisAngle = 50 #in deg
SemiLatusRectum = 60 # in mm,SemiLatusRectum =  2*focal length of the "mother parabola"
Mirror5 = mmirror.MirrorParabolic(SemiLatusRectum, offAxisAngle, msupp.SupportRound(25) ) 

OpticsList = [Mirror1, Mirror2, Mirror3, Mirror4, Mirror5]

OpticalChain = mp.OEPlacement(Divergence, SourceSize, NumberRays, OpticsList, DistanceList, AngleList, IncidencePlaneAngleList)


#%% detector
ReflectionNumber = -1 # place detector after mirror number, -1 being the last

ManualDetector = False # just have it automatically placed perpendicular to the final ray bundle at distance DistanceMirrorDetector from last mirror
DistanceDetector = SemiLatusRectum/(1+np.cos(offAxisAngle/180*np.pi))
AutoDetectorDistance = False  # have ART find the optimal detector distance
            


#%%#####################################################################################################################
""" Do THE ACTUAL RAYTRACING calculations """
""" You can loop over this block to simulate a range of optical setups, e.g. a varying distance """

exec( open("ART/DoART.py").read() )

#%%
"""  THE built-in plot options """

exec( open("ART/DoPLOTS.py").read() )

########################################################################################################################        
   





