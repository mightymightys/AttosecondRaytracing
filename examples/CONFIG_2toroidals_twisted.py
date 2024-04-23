# -*- coding: utf-8 -*-
"""
Created in Jan 2023

@author: Haessler
"""
#%% Modules
#import copy
import numpy as np
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleMask as mmask
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleProcessing as mp
from ART.ARTmain import main

#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 50e-3/2, # half-angle in rad, 0 for a plane wave
    'SourceSize' : 0,       # diameter in mm, 0 for a point source
    'Wavelength' : 50e-6,   # 50 nm
    'DeltaFT'    : 0.5,     # in fs
    'NumberRays' : 1000     # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
Description = '2 toroidal mirrors in f-D-f config, i.e. approx. collimation, propagation, and the refocus, with the 2nd toroidal rotating the plane of incidence.'

# a mask that sits 500mm from the source and limits the incident ray bundle to a 41mrad divergence full-angle
Mask = mmask.Mask(msupp.SupportRoundHole(30,41e-3/2*500,0,0))

#2 equal toroidal mirrors
Support = msupp.SupportRectangle(200,30)
FocalLength = 600
ToroidalIncidence = 80 #in deg
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(FocalLength, ToroidalIncidence)
ToroidalMirror = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)

# preparing the optical chain
OpticsList = [Mask, ToroidalMirror, ToroidalMirror ]

IncidenceAngleList = [0, ToroidalIncidence, -ToroidalIncidence] 
DistanceList = [500, FocalLength-500, FocalLength]
# loop over the incidence-plane rotation by the 2nd toroidal mirror:
IncidencePlaneAngleList = [0, 0, np.linspace(-90, 90, 10) ]

# produce a png-image of each of the varied optical chains ?
render = False
OpticalChainList =  mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description)

#OpticalChainList[0].quickshow()



#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1,       # analyse ray bundle after this optical element (index of you optics-list above)
    'ManualDetector' : False,      # do not place the detector manually, i.e. let ART do it automatically at the distance 'DistanceDetector' and perpendicular to the central ray
    'DistanceDetector' : FocalLength ,     # set the detector at this distance from the selected optical element
    'AutoDetectorDistance' : True, # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"         # metric for which to optimize the detector position
}

#%%
"""
ANALYSIS OPTIONS
You can delete lines that don't interest you out of this disctionary.
They will be automatically completed by default values.
All plot_*-options are False by default.
"""
AnalysisOptions = {
    'verbose': True,         # print intermediate results and info in the console?

    'plot_Render': True,    # render optical elements and rays?
    
    'DrawAiryAndFourier': True,  # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,          # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': True,     # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': False,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
    'plot_IncidenceGraph': False,    # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?

    'plot_DelayMirrorProjection': False,      # produce a plot of the ray delays at the detector projected onto the surface of the preceding mirror?
    'plot_IntensityMirrorProjection': False,  # produce a plot of the ray intensities at the detector projected onto the surface of the preceding mirror?
    'plot_IncidenceMirrorProjection': False,  # produce a plot of the ray incidence angles at the detector projected onto the surface of the preceding mirror?

    'save_results': True  #save the simulation results to disk, to analyse later
}

#%%
"""
TO RUN A SIMULATION, TYPE IN AN ANACONDA-PROMPT (or equivalent):
    python ARTmain.py <name_of_this_configuration_file>
    
    TO LOAD THE SAVED kept_data DICTIONARY FROM DISK AND USE IT FURTHER, DO
    import ART.ModuleProcessing as mp
    kept_data = mp.load_compressed(filename)
    
    
IF PREFER RUNNING THIS CONFIG-SCRIPT DIRECTLY, call the main function of the
ARTmain.py-program from here like so:
"""
if __name__ == "__main__":
    kept_data = main(OpticalChainList,
                    SourceProperties,
                    DetectorOptions,
                    AnalysisOptions)

# """ IN EITHER CASE, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO: """
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD (fs)")
