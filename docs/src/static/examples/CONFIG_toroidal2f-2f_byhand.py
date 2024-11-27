# -*- coding: utf-8 -*-
"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules
#import copy

import numpy as np
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleMask as mmask
from ART.ARTmain import main
import ARTcore.ModuleOpticalElement as moe
import ARTcore.ModuleSource as msource
import ARTcore.ModuleOpticalChain as moc
#import ART.ModuleProcessing as mp


#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 15e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 50e-6, # 50 nm
    'DeltaFT'    : 0.5,   # in fs
    'NumberRays' : 1000   # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
description = 'single toroidal in 2f-2f config,\n set up "byhand", with the mirror sitting at the origin and the source at [Sx,Sy,Sz]'

#% the toroidal mirror
DimX = 120
DimY = 30
Support = msupp.SupportRectangle(DimX,DimY)
    
AngleIncidence = 80
Focal = 300
R, r = mmirror.ReturnOptimalToroidalRadii(Focal, AngleIncidence)
Mirror = mmirror.MirrorToroidal(R, r, Support)

# set up the optical element corresponding to the mirror:
# it sits at the origin, the normal and majoraxis point in the z and x direction, respectively
Element1 = moe.OpticalElement(Mirror, np.array([0,0,0]), np.array([0,0,1]), np.array([1,0,0]))

MisalignAngle = 0.00 #deg

#The source sits at a point 2*Focal away at an angle AngleIncidence with the z-axis
Sx = 2*Focal *np.sqrt(1-np.cos(np.deg2rad(AngleIncidence+MisalignAngle))**2)
Sy = 0
Sz = 2*Focal * np.cos(np.deg2rad(AngleIncidence+MisalignAngle))
SourcePoint = np.array([Sx,Sy,Sz])
SourceRayList = msource.PointSource(SourcePoint, -SourcePoint, SourceProperties["Divergence"], SourceProperties["NumberRays"])
SourceRayList = msource.ApplyGaussianIntensityToRayList(SourceRayList, 1/np.e**2)

#% creating the optical chain
OpticalChainList = moc.OpticalChain(SourceRayList, [Element1], description, loop_variable_name = "incidence misalignment (deg)", loop_variable_value =MisalignAngle)
    

#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1, # analyse ray bundle after last optical element
    'ManualDetector' : False, # do not place the detector manually, i.e. let ART do it automatically
    'DistanceDetector' : 2*Focal ,  # set the detector at this distance from the last optical element
    'AutoDetectorDistance' : False,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"   # metric for which to optimize the detector position
}

#%%
""" Analysis options """

AnalysisOptions = {
    'verbose': False,           # print intermediate results and info in the console?

    'plot_Render': True,         # render optical elements and rays?

    'DrawAiryAndFourier': True,    # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,          # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': False,  # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': True,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
    'plot_IncidenceGraph': False,    # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?

    'plot_DelayMirrorProjection': False,      # produce a plot of the ray delays at the detector projected onto the surface of the preceding mirror?
    'plot_IntensityMirrorProjection': False, # produce a plot of the ray intensities at the detector projected onto the surface of the preceding mirror?
    'plot_IncidenceMirrorProjection': False,  # produce a plot of the ray incidence angles at the detector projected onto the surface of the preceding mirror?

    'save_results': False        #save the simulation results to disk, to analyse later
}

#%%
"""
TO RUN A SIMULATION, TYPE IN AN ANACONDA-PROMPT (or equivalent):
    python ARTmain.py <name_of_this_configuration_file>
    
    TO LOAD THE SAVED kept_data DICTIONARY FROM DISK AND USE IT FURTHER, DO
    import ART.ModuleProcessing as mp
    kept_data = mp.load_compressed(filename)
    
    
IF WANT TO RUN THIS CONFIG-SCRIPT DIRECTLY, call the main function of the ARTmain.py-program from here:

from ARTmain import main
kept_data = main(OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions)
 
  
#    THEN, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO:
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
ax.set_ylabel("DurationSD")
"""

if __name__ == "__main__":
    kept_data = main(OpticalChainList,
                    SourceProperties,
                    DetectorOptions,
                    AnalysisOptions)
