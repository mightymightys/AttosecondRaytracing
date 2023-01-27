"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules

import numpy as np

import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp


#%%########################################################################
""" User cell """
verbose = True  	# print intermediate results and info in the console?

FindBestSpot = False      	# look for detector distance with smallest spotsize?
FindBestDuration = False  	# look for detector distance with shortest duration?
FindOptDistance  = False  	# look for optimal detector distance with "highest intensity"?

RayGraph = True  	# render optical elements and rays, and how many rays to render? 
maxRaysToRender = 200

DrawAiryAndFourier = True     # Draw Airy spot and Fourier-limited duration in the following plots?
SpotDiagram = False           # produce an interactive spot diagram without color coding the spots?	
DelaySpotDiagram = True       # produce an interactive spot diagram with ray delays color coded?
IntensitySpotDiagram = False  # produce an interactive spot diagram with ray intensities color coded?
IncidenceSpotDiagram = False  # produce an interactive spot diagram with ray incidence angles color coded?
	
DelayGraph = False      # produce an interactive spot diagram with delays in 3rd dimension?
IntensityGraph = True   # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
IncidenceGraph = False  # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?
	
DelayMirrorProjection = True      # produce a plot of the ray delays at the detector projected onto the mirror surface?
IntensityMirrorProjection = False # produce a plot of the ray intensities at the detector projected onto the mirror surface?
IncidenceMirrorProjection = False # produce a plot of the ray incidence angles at the detector projected onto the mirror surface?


#%%
""" try a cylindrical mirror """

#% source properties
Divergence = 15e-3  # half-angle in rad; 
SourceSize = 0e-3   # in mm
Wavelength = 50e-6  # 50 nm, only plays a role for the Airy spot size show in some plots as a reference
DeltaFT = 0.5       # 500 as, only appears as a reference in some plots
NumberRays = 1000

#% Mirror
Support = msupp.SupportRectangle(60,30)

Focal = 100
Mirror = mmirror.MirrorCylindrical( 2*Focal, Support)

#% creating the optical chain
OpticsList = [Mirror]
IncidenceAngleList = [10] #in deg
IncidencePlaneAngleList = [0]
DistanceList = [2*Focal]
    
OpticalChain = mp.OEPlacement(Divergence, SourceSize, NumberRays, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList)
  

#% detector
ReflectionNumber = -1 # place detector after mirror number, -1 being the last

ManualDetector = False # just have it automatically placed perpendicular to the final ray bundle at distance DistanceMirrorDetector from last mirror
DistanceDetector = 2*Focal
AutoDetectorDistance = False  # have ART find the optimal detector distance

    
#%%#####################################################################################################################
""" Do THE ACTUAL RAYTRACING calculations """
""" You can loop over this block to simulate a range of optical setups, e.g. a varying distance """

exec( open("ART/DoART.py").read() )

#%%   
"""  Go through the BUILT-IN PLOT OPTIONS """

exec( open("ART/DoPLOTS.py").read() )
########################################################################################################################
        
   



