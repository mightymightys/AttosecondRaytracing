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
    'Divergence' : 0, # half-angle in rad, 0 for a plane wave
    'SourceSize' : 50,      # diameter in mm, 0 for a point source
    'Wavelength' : 800e-6,  # 800 nm
    'DeltaFT'    : 2.7,     # in fs
    'NumberRays' : 1000     # launch 1000 rays in the beginning 
}


#%%
""" OPTICAL SETUP """
Description = "A 90° off-axis parabola with a hole, illuminated by a plane wave."

# Mirror
Support =  msupp.SupportRoundHole(30,5,10,5)
#Support =  msupp.SupportRound(50)

offAxisAngle = 90 #in deg
FocalEffective = 100 # in mm
Parabola = mmirror.MirrorParabolic(FocalEffective, offAxisAngle, Support) 

# creating the optical chain
OpticsList = [Parabola]
DistanceList = [200] 
IncidenceAngleList = [0.00] # in degrees, for the parabola this angle is relative to the mother parabola symmetry axis


OpticalChainList = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, Description = Description)

"""
Now that the optical chain has been automatically created with perfect alignment, you can mess with it:
As a demonstration, mis-align the parabola slightly out of the incidence plane, i.e. rotate its normal about its major axis,
which is what the rotate_roll_by(angle in degress)-method of the optical elements does:
"""
ParabolaOE = OpticalChainList.optical_elements[0] #pick out the parabola (it's the only element in the list anyway)
turnby = np.rad2deg(50e-6) #in deg
ParabolaOE.rotate_roll_by(turnby)  

# This is enough, we have modified directly a property of the OpticalChain-object. 



#%%
""" detector parameters """

DetectorOptions = {
    'ReflectionNumber' : -1,       # analyse ray bundle after this optical element (index of you optics-list above)
    'ManualDetector' : False,      # do not place the detector manually, i.e. let ART do it automatically at the distance 'DistanceDetector' and perpendicular to the central ray
    'DistanceDetector' : FocalEffective ,     # set the detector at this distance from the selected optical element
    'AutoDetectorDistance' : False, # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"          # metric for which to optimize the detector position
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

    'plot_Render': True,    # render optical elements and rays

    'DrawAiryAndFourier': True,  # Draw Airy spot and Fourier-limited duration in the following plots?
    
    'plot_SpotDiagram': False,          # produce an interactive spot diagram without color coding the spots?
    'plot_DelaySpotDiagram': True,     # produce an interactive spot diagram with ray delays color coded?
    'plot_IntensitySpotDiagram': False, # produce an interactive spot diagram with ray intensities color coded?
    'plot_IncidenceSpotDiagram': False, # produce an interactive spot diagram with ray incidence angles color coded?

    'plot_DelayGraph': False,        # produce an interactive spot diagram with delays in 3rd dimension?
    'plot_IntensityGraph': True,     # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
    'plot_IncidenceGraph': False,    # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?

    'plot_DelayMirrorProjection': False,      # produce a plot of the ray delays at the detector projected onto the surface of the preceding mirror?
    'plot_IntensityMirrorProjection': False,  # produce a plot of the ray intensities at the detector projected onto the surface of the preceding mirror?
    'plot_IncidenceMirrorProjection': False,  # produce a plot of the ray incidence angles at the detector projected onto the surface of the preceding mirror?

    'save_results': False  #save the simulation results to disk, to analyse later
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
# kept_data = main(OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions)
    

# """ IN EITHER CASE, YOU CAN PLOT RESULTS AGAINST e.g. THE LOOP-VARIABLE LIKE SO: """
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.scatter([x.loop_variable_value for x in kept_data["OpticalChain"]], kept_data["DurationSD"])
# ax.set_xlabel(kept_data["OpticalChain"][0].loop_variable_name)
# ax.set_ylabel("DurationSD (fs)")
if __name__ == "__main__":
    kept_data = main(OpticalChainList,
                    SourceProperties,
                    DetectorOptions,
                    AnalysisOptions)
