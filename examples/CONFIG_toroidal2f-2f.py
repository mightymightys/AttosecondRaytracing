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
from ART.ARTmain import main

#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 30e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 50e-6, # 50 nm
    'DeltaFT'    : 0.5,   # in fs
    'NumberRays' : 1000   # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
""" first create a perfectly aligned optical chain """
Description = "single toroidal or ellipsoidal in 2f-2f config, possibly misaligned "

Support = msupp.SupportRectangle(300,50)

Focal = 500
MirrorIncidence = 80 #in deg
# Toroidal
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal,MirrorIncidence)
Mirror = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)
## or Ellipsoidal
# SemiMajorAxis, SemiMinorAxis = mmirror.ReturnOptimalEllipsoidalAxes(2*Focal, MirrorIncidence)
# Mirror = mmirror.MirrorEllipsoidal(SemiMajorAxis, SemiMinorAxis, Support)

# create optical chain
OpticsList = [Mirror]
DistanceList = [2*Focal]
IncidenceAngleList = [MirrorIncidence] #in deg
AlignedOpticalChain = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, Description = Description)

""" then make a list of optical chains, each one with some modification with respect to the aligned one """
OEindx = 0
loop_axis = "roll"
loop_variable_values = np.linspace(-0.5, 0.5, 11)
OpticalChainList = AlignedOpticalChain.get_OE_loop_list(OEindx, loop_axis, loop_variable_values )



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

    'plot_Render': False,         # render optical elements and rays?

    'DrawAiryAndFourier': True,    # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,       # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': False,  # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': False,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
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
"""
# from ARTmain import main
# kept_data = main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions)
    
# #%%    THEN, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO:
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD")

# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], [1e3*x for x in kept_data["SpotSizeSD"]])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("SpotSizeSD ($\mu m$)")

if __name__ == "__main__":
    kept_data = main(OpticalChainList,
                    SourceProperties,
                    DetectorOptions,
                    AnalysisOptions)
