"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules

import numpy as np

import ART.ModuleMirror as mmirror
import ART.ModuleMask as mmask
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp
import ART.ModuleGeometry as mgeo

#%%########################################################################
""" User cell """
verbose = True  	# print intermediate results and info in the console?

FindBestSpot = False      	# look for detector distance with smallest spotsize?
FindBestDuration = False  	# look for detector distance with shortest duration?
FindOptDistance  = False  	# look for optimal detector distance with "highest intensity"?

RayGraph = True     	    # render optical elements and rays, and how many rays to render? 
maxRaysToRender = 200

DrawAiryAndFourier = True     # Draw Airy spot and Fourier-limited duration in the following plots?
SpotDiagram = False           # produce an interactive spot diagram without color coding the spots?	
DelaySpotDiagram = True       # produce an interactive spot diagram with ray delays color coded?
IntensitySpotDiagram = False  # produce an interactive spot diagram with ray intensities color coded?
IncidenceSpotDiagram = False  # produce an interactive spot diagram with ray incidence angles color coded?
	
DelayGraph = False      # produce an interactive spot diagram with delays in 3rd dimension?
IntensityGraph = False   # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
IncidenceGraph = False  # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?
	
DelayMirrorProjection = False     # produce a plot of the ray delays at the detector projected onto the mirror surface?
IntensityMirrorProjection = False # produce a plot of the ray intensities at the detector projected onto the mirror surface?
IncidenceMirrorProjection = False # produce a plot of the ray incidence angles at the detector projected onto the mirror surface?

#%%

"""single parabola"""
#% source properties
Divergence = 0      # half-angle in rad, 0 creates a plane wave
SourceSize = 60     # in mm
Wavelength = 800e-6 # 800nm in mm 
DeltaFT = 3         # Fourier-limited pulse duration in fs
NumberRays = 1000

#% mirror
Support =  msupp.SupportRound(40)

offAxisAngle = 5 #in deg
SemiLatusRectum = 3000 # in mm,SemiLatusRectum =  2*focal length of the "mother parabola"
Parabola1 = mmirror.MirrorParabolic(SemiLatusRectum, offAxisAngle, Support) 
Parabola1bis = mmirror.MirrorParabolic(SemiLatusRectum, -offAxisAngle, Support)
FocalEffective1 = SemiLatusRectum/(1+np.cos(offAxisAngle/180*np.pi))

#Spherical = mmirror.MirrorSpherical(3000, Support) 
#FocalEffective1 = 1500

offAxisAngle = 90 #in deg
SemiLatusRectum = 100 # in mm
Parabola2 = mmirror.MirrorParabolic(SemiLatusRectum, offAxisAngle, Support) 
FocalEffective2 = SemiLatusRectum/(1+np.cos(offAxisAngle/180*np.pi))

PM = mmirror.MirrorPlane(Support) 

#% creating the optical chain
OpticsList = [Parabola1, PM, Parabola1bis,  Parabola2]
DistanceList = [FocalEffective1, FocalEffective1, FocalEffective1, 2000] 
IncidenceAngleList = [0, 5, -5, 0] # in degrees, for the parabola this angle is relative to the mother parabola symmetry axis

#OpticsList = [Spherical, PM, Spherical,  Parabola2]
#DistanceList = [FocalEffective1, FocalEffective1, FocalEffective1, 10000] 
#IncidenceAngleList = [-2.5, 5, 2.5, 0] # in degrees, for the parabola this angle is relative to the mother parabola symmetry axis


OpticalChain = mp.OEPlacement(Divergence, SourceSize, NumberRays, OpticsList, DistanceList, IncidenceAngleList)

"""
Now that the optical chain has been automatically created with perfect alignment, we can mess with it:
As a demonstration, mis-align the parabola slightly out of the incidence plane, i.e. rotate its normal about its major axis:
"""
#ParabolaOE = OpticalChain[1]
#turnby = 0.00 # deg
#ParabolaOE.normal = mgeo.RotationAroundAxis(ParabolaOE.majoraxis, turnby/180*np.pi, ParabolaOE.normal)    


#% detector
ReflectionNumber = -4 # analyse ray bundle after last optical element

ManualDetector = False # do not place the detector manually, i.e. let ART do it automatically
DistanceDetector = FocalEffective2  # set the detector at this distance from the last optical element
AutoDetectorDistance = True  # search for the optimal detector distance and shift the detector there ?


#%%#####################################################################################################################
""" Do THE ACTUAL RAYTRACING calculations """
""" You can loop over this block to simulate a range of optical setups, e.g. a varying distance """

exec( open("ART/DoART.py").read() )


#%%   
"""  Go through the BUILT-IN PLOT OPTIONS """

exec( open("ART/DoPLOTS.py").read() )
########################################################################################################################