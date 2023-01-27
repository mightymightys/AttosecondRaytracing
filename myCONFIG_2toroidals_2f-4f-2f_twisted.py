"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules
import numpy as np
import pandas as pd
#from tqdm import tqdm

import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp


#%%########################################################################
""" User cell """
verbose = True  	# print intermediate results and info in the console?

FindBestSpot = False      	# look for detector distance with smallest spotsize?
FindBestDuration = False  	# look for detector distance with shortest duration?
FindOptDistance  = False  	# look for optimal detector distance with "highest intensity"?

RayGraph = False       	# render optical elements and rays, and how many rays to render? 
maxRaysToRender = 200

DrawAiryAndFourier = True     # Draw Airy spot and Fourier-limited duration in the following plots?
SpotDiagram = False           # produce an interactive spot diagram without color coding the spots?	
DelaySpotDiagram = True       # produce an interactive spot diagram with ray delays color coded?
IntensitySpotDiagram = False  # produce an interactive spot diagram with ray intensities color coded?
IncidenceSpotDiagram = False  # produce an interactive spot diagram with ray incidence angles color coded?
	
DelayGraph = False      # produce an interactive spot diagram with delays in 3rd dimension?
IntensityGraph = False   # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
IncidenceGraph = False  # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?
	
DelayMirrorProjection = True      # produce a plot of the ray delays at the detector projected onto the mirror surface?
IntensityMirrorProjection = False # produce a plot of the ray intensities at the detector projected onto the mirror surface?
IncidenceMirrorProjection = False # produce a plot of the ray incidence angles at the detector projected onto the mirror surface?


#%%
"""double twisted toroidal config """

#% source properties
Divergence = 30e-3/2  # half-angle in rad
SourceSize = 2e-3    # in mm
Wavelength = 40e-6  # 50 nm, only Airy spot size show in some plots
DeltaFT = 0.5       # 500 as, only appears as a reference in some plots
NumberRays = 2000   # launch 1000 rays in the beginning 

#% mirrors
Support = msupp.SupportRectangle(120,30)
Focal = 300
ToroidalIncidence = 80 #in deg
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal, ToroidalIncidence)
ToroidalMirror = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius, Support)

#PlaneMirror = mmirror.MirrorPlane(Support)

#% creating the optical chain
OpticsList = [ToroidalMirror, ToroidalMirror]
misalignby = 0.00 #deg
IncidenceAngleList = [ToroidalIncidence+misalignby, -ToroidalIncidence-misalignby] #in deg

DistanceList = [2*Focal, 4*Focal] # 2f-4f-2f config, i.e. with intermediate focus
#DistanceList = [Focal, Focal] # f-f-f config, i.e. with intermediate "collimated beam" and optimal mirror distance

#% detector
ReflectionNumber = -1        # analyse ray bundle after last optical element
ManualDetector = False       # do not place the detector manually, i.e. let ART do it automatically
DistanceDetector = 2*Focal  # set the detector at this distance from the last optical element
AutoDetectorDistance = True  # then search for the optimal detector distance and shift the detector there?


#%%
#%%
# initialize lists to contains all the results we want to keep from the loop iterations
RotationS = [] 
RayListHistoryS = []
OpticalChainS = []
DetectorS = []
OptDistanceS = []
OptSizeSpotS = []
OptDurationS = []

#for Rotation in tqdm(range(0, 350, 10)):
#for Rotation in tqdm(range(0, 185, 5)):
for Rotation in [180]:
    
    if Rotation >= 180: 
        RotationBis = Rotation-180
        IncidenceAngleList = [ToroidalIncidence+misalignby, ToroidalIncidence+misalignby]
    else:
        RotationBis = Rotation
        
    IncidencePlaneAngleList = [0,RotationBis]
    
    OpticalChain = mp.OEPlacement(Divergence, SourceSize, NumberRays, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList)
    
    #%%#####################################################################################################################
    """ Do THE ACTUAL RAYTRACING calculations """
    """ You can loop over this block to simulate a range of optical setups, e.g. a varying distance """
    ########################################################################################################################
    
    exec( open("ART/DoART.py").read() )

    # append the lists to save the results of all iterations
    RotationS.append(Rotation)
    RayListHistoryS.append(RayListHistory)
    OpticalChainS.append(OpticalChain)
    DetectorS.append(Detector)
    OptDistanceS.append(OptDistance)
    OptSizeSpotS.append(OptSizeSpot)
    OptDurationS.append(OptDuration)

   
    # mlab.options.offscreen = True
    # rayRenderFig = mplots.RayRenderGraph(RayListHistory,OpticalChain,Detector.get_distance()*1.5 ,maxRaysToRender)
    # name = "2f-4f-2f_for_{foo:03d}deg.png"
    # mlab.savefig(name.format(foo=int(Rotation)), figure=rayRenderFig, magnification=3)


#%% create a pandas DataFrame to contain all the saved lists with descriptive column titles, which can then be saved to an efficient hdf file
resultsDF = pd.DataFrame({'Mirror Rotation [deg]':RotationS, 'RayListHistory':RayListHistoryS, 'OpticalChain':OpticalChainS, 'OptDetector':DetectorS, 'OptDetDistance [mm]':OptDistanceS, 'OptSpotSize [micron]':[x*1000 for x in OptSizeSpotS], 'OptDuration [fs]':OptDurationS})
scanfig = resultsDF.plot(kind='scatter',x='Mirror Rotation [deg]',y='OptDuration [fs]')


#%%
"""possibly store results"""
# name = "2f-4f-2f_for_"+str(int(Divergence*1000))+"mrad.h5"
# resultsDF.to_hdf(name, 'resultsDF'+str(int(Divergence*1000))+'mrad')


#%%#####################################################################################################################
# Select the setup that gives the shortest duration at its best focus    
MinDurationRow = resultsDF['OptDuration [fs]'].idxmin()

OpticalChain = resultsDF.at[MinDurationRow,'OpticalChain']
RayListHistory = resultsDF.at[MinDurationRow,'RayListHistory']
RayListAnalysed = RayListHistory[ReflectionNumber]
Detector = resultsDF.at[MinDurationRow,'OptDetector']

#%
"""  Go through the BUILT-IN PLOT OPTIONS """

exec( open("ART/DoPLOTS.py").read() )
########################################################################################################################
        
   



