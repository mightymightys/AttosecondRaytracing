# -*- coding: utf-8 -*-
"""
Created in Jan 2023

@author: Haessler
"""
#%% Modules
#import copy
import numpy as np
import ART.ModuleMirror as mmirror
#import ART.ModuleMask as mmask
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp


#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 0.8e-3,  # half-angle in rad, 0 for a plane wave
    #'Divergence' : 3e-3, 
    #'SourceSize' : 0,
    'SourceSize' : 120e-3,   # diameter in mm, 0 for a point source
    #'SourceSize' : 30e-3,
    'Wavelength' : 50e-6,   # 40 nm
    'DeltaFT'    : 0.5,     # in fs
    'NumberRays' : 1000     # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
Description = "3 three toroidal mirrors to microfocus, as in [Polettoet al., Opt. Express 21, 13040 (2013)]:\n\
    f_1-d-f_2+2*f_3-2*f_3 config, i.e. approx. collimation, propagation, 'micro-refocus', \n\
    and then approx. 2f-2f reimaging of that small focus."

# 3 toroidal mirrors with same support
Support = msupp.SupportRectangle(100,30)

AngleIncidence = 80 #in deg

Focal1 = 1500
#Focal1 = 300
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal1, AngleIncidence)
ToroidalMirror1 = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)

Focal2 = 150
#Focal2 = 100
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal2, AngleIncidence)
ToroidalMirror2 = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)

Focal3 = 314.5
#Focal3 = 150
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal3, AngleIncidence)
ToroidalMirror3 = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)

# creating the optical chain
OpticsList = [ToroidalMirror1, ToroidalMirror2, ToroidalMirror3]
IncidenceAngleList = [AngleIncidence,-AngleIncidence, AngleIncidence] #in deg
IncidencePlaneAngleList = [0,0,0]

# now make a list of optical chains with varying distance between the first and second mirror, i.e. the quasi-collimated section
OpticalChainList = []
loop_variable_name = "length of collimated section (mm)"
loop_variable =  np.linspace(500, 30000, 11)

for distance in loop_variable:
    DistanceList = [Focal1, distance, Focal2+660]
    #DistanceList = [Focal1, Distance, Focal2+321]

    ModifiedOpticalChain = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description = Description)
    ModifiedOpticalChain.loop_variable_name = loop_variable_name
    ModifiedOpticalChain.loop_variable_value = distance

    OpticalChainList.append(ModifiedOpticalChain)

"""
At the end of this OPTICAL SETUP block, we have a list of optical-chain-objects,
called OpticalChainList. Good.
"""


#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1,       # analyse ray bundle after this optical element (index of you optics-list above)
    'ManualDetector' : False,      # do not place the detector manually, i.e. let ART do it automatically at the distance 'DistanceDetector' and perpendicular to the central ray
    'DistanceDetector' : 600,     # set the detector at this distance from the selected optical element
    #'DistanceDetector' : 282,
    'AutoDetectorDistance' : True, # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"         # metric for which to optimize the detector position
}

#%%
"""
ANALYSIS OPTIONS
"""
AnalysisOptions = {
    'verbose': False,         # print intermediate results and info in the console?

    'plot_Render': False,    # render optical elements and rays, and how many rays to render (more than 200 gets very slow)?
    'maxRaysToRender': 150,

    'DrawAiryAndFourier': True,  # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,          # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': False,     # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': True,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
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
# from ARTmain import main
# kept_data = main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions)
    
# """ THEN, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO: """
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD")
