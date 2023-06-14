"""
Created in October 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules

import numpy as np
import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp
import ART.ModuleMask as mmask


#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 50e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 80e-6, # 40 nm
    'DeltaFT'    : 0.5,   # in fs
    'NumberRays' : 1000   # launch 1000 rays in the beginning 
}

#%%
""" OPTICAL SETUP """
""" first create a perfectly aligned optical chain as usual """
Description = "2 toroidal mirrors in f-d-f config, i.e. approx. collimation, propagation, and the refocus "

#% Mirrors
Support = msupp.SupportRectangle(150, 32)

Focal = 500
AngleIncidence = 80 #in deg
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal, AngleIncidence)
ToroidalMirror = mmirror.MirrorToroidal(OptimalMajorRadius, OptimalMinorRadius, Support)

# a mask
#placed 400 after source, this selects 35mrad: 
SupportMask = msupp.SupportRoundHole(Radius=20, RadiusHole=14/2, CenterHoleX=0, CenterHoleY=0) 
Mask = mmask.Mask(SupportMask)

# create the optical chain
OpticsList = [Mask, ToroidalMirror, ToroidalMirror]
IncidenceAngleList = [0, AngleIncidence, -AngleIncidence] #in deg
IncidencePlaneAngleList = [0, 0, 0]
# loop over the distance between the 2 toroidal mirrors:
DistanceList = [400, Focal-400, np.linspace(Focal-200, Focal+200, 11)] #one element is an array of values, so a list of optical chains will be created

# produce a png-image of each of the varied optical chains ?
render = False
OpticalChainList =  mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList, Description, render)

#OpticalChainList[0].quickshow()


#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1, # analyse ray bundle after last optical element
    'ManualDetector' : False, # do not place the detector manually, i.e. let ART do it automatically
    'DistanceDetector' : Focal ,  # set the detector at this distance from the last optical element
    'AutoDetectorDistance' : True,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"   # metric for which to optimize the detector position
}

#%%
""" Analysis options """

AnalysisOptions = {
    'verbose': True,           # print intermediate results and info in the console?

    'plot_Render': True,         # render optical elements and rays, and how many rays to render (more than 200 gets very slow)?
    'maxRaysToRender': 150,

    'DrawAiryAndFourier': True,    # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,          # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': True,  # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': False,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
    'plot_IncidenceGraph': False,    # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?

    'plot_DelayMirrorProjection': False,      # produce a plot of the ray delays at the detector projected onto the surface of the preceding mirror?
    'plot_IntensityMirrorProjection': False, # produce a plot of the ray intensities at the detector projected onto the surface of the preceding mirror?
    'plot_IncidenceMirrorProjection': False,  # produce a plot of the ray incidence angles at the detector projected onto the surface of the preceding mirror?

    'save_results': True        #save the simulation results to disk, to analyse later
}

#%%
"""
TO RUN A SIMULATION, TYPE IN AN ANACONDA-PROMPT (or equivalent):
    python ARTmain.py <name_of_this_configuration_file>
    
    TO LOAD THE SAVED kept_data DICTIONARY FROM DISK AND USE IT FURTHER, DO
    import ART.ModuleProcessing as mp
    kept_data = mp.load_compressed(filename)
    
    
IF WANT TO RUN THIS CONFIG-SCRIPT DIRECTLY, call the main function of the ARTmain.py-program from here:
"""
# from ARTmain import main
# kept_data = main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions)
    
# #    THEN, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO:
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD")



from ARTmain import main
kept_data = main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions)
