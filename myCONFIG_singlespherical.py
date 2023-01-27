"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules
#import numpy as np
import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp

#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 15e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 40e-6, # 780 nm, typical central WL for postcompressed Ti:Sa laser; only plays a role for the Airy spot size show in some plots as a reference
    'DeltaFT'    : 0.5,    # half-period of 780-nm wave (here we care for spatio-temporal distortions on the optical-cycle rather than envelope scale)
    'NumberRays' : 1000    # launch 1000 rays in the beginning 
}



#%%
""" OPTICAL SETUP """    
Description = " Single spherical mirror, near normal incidence, MPQ-like """
#% mirror
Focal = 100
SphMirror = mmirror.MirrorSpherical(2*Focal, msupp.SupportRound(50))

#% creating the optical chain
OpticsList = [SphMirror]

DistanceList = [1000]  # in mm
IncidenceAngleList = [2]  # in degrees

OpticalChain = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, Description = Description)

#OpticalChain.quickshow()


#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1, # analyse ray bundle after last optical element
    'ManualDetector' : False, # do not place the detector manually, i.e. let ART do it automatically
    'DistanceDetector' : Focal,  # set the detector at this distance from the last optical element
    'AutoDetectorDistance' : True,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
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
    'plot_DelaySpotDiagram': True,  # produce an interactive spot diagram with ray delays color coded?
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
    
IF WANT TO RUN THIS CONFIG-SCRIPT DIRECTLY, call the main function of the ARTmain.py-program from here:
    from ARTmain import main
    kept_data = main(OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions, save_file_name)
"""