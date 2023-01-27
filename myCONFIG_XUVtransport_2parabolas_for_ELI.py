# -*- coding: utf-8 -*-
"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules
#import copy
import numpy as np
import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp


#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 50e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 40e-6, # 40 nm
    'DeltaFT'    : 0.5,   # in fs
    'NumberRays' : 1000   # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
""" first create a perfectly aligned optical chain as usual """
Description = "2 equal large-off-axis-angle parabolas for collimation and refocusing, with plane mirror in the middle "

# parabolas
Support = msupp.SupportRound(25)

offAxisAngle = 140 #in deg
FocalEffective = 500 # in mm
SemiLatusRectum = FocalEffective*(1+np.cos(offAxisAngle/180*np.pi)) # in mm
Parabola = mmirror.MirrorParabolic(SemiLatusRectum, offAxisAngle, Support) 
ParabolaBis = mmirror.MirrorParabolic(SemiLatusRectum, -offAxisAngle, Support) 

# a plane mirror
PM = mmirror.MirrorPlane(msupp.SupportRectangle(100,50)) 

# create the optical chain
OpticsList = [Parabola, PM, ParabolaBis]
DistanceList = [FocalEffective, 500, 500] 
IncidenceAngleList = [offAxisAngle, -70, 0] # in deg, for the parabola this angle is relative to the mother parabola symmetry axis

AlignedOpticalChain = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, Description = Description)

# """ the make a list of optical chains, each one with some modification with respect to the aligned one """
# OpticalChainList = []
# loop_variable_name = "roll misalignement angle (deg)"
# loop_variable =  np.linspace(-0.5, 0.5, 11)

# for angle in loop_variable:
#     # always start with a fresh deep-copy the AlignedOpticalChain, to then modify it and append it to the list
#     ModifiedOpticalChain = AlignedOpticalChain.copy_chain()
#     ModifiedOpticalChain.loop_variable_name = loop_variable_name
#     ModifiedOpticalChain.loop_variable_value = angle
    
#     #pick out the optical element you want to mess with and modify its alignment
#     ModifiedMirrorOE = ModifiedOpticalChain.optical_elements[0] 
#     ModifiedMirrorOE.rotate_roll_by(angle)
    
#     #append the modified optical chain to the list
#     OpticalChainList.append(ModifiedOpticalChain)

OpticalChainList = [AlignedOpticalChain]

#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1, # analyse ray bundle after last optical element
    'ManualDetector' : False, # do not place the detector manually, i.e. let ART do it automatically
    'DistanceDetector' : FocalEffective ,  # set the detector at this distance from the last optical element
    'AutoDetectorDistance' : False,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"   # metric for which to optimize the detector position
}

#%%
""" Analysis options """

AnalysisOptions = {
    'verbose': True,           # print intermediate results and info in the console?

    'plot_Render': False,         # render optical elements and rays, and how many rays to render (more than 200 gets very slow)?
    'maxRaysToRender': 150,

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
"""
# from ARTmain import main
# kept_data = main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions)
    
# #    THEN, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO:
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD")
