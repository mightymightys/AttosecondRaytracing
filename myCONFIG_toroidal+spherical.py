"""
Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
#%% Modules

import numpy as np
#import addcopyfighandler  #if you want to be able to copy figure contents to the clipboard by ctrl-c, Only for windows though.
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

RayGraph =  False	    # render optical elements and rays, and how many rays to render? 
maxRaysToRender = 100

DrawAiryAndFourier = True     # Draw Airy spot and Fourier-limited duration in the following plots?
SpotDiagram = False           # produce an interactive spot diagram without color coding the spots?	
DelaySpotDiagram = True       # produce an interactive spot diagram with ray delays color coded?
IntensitySpotDiagram = False  # produce an interactive spot diagram with ray intensities color coded?
IncidenceSpotDiagram = False  # produce an interactive spot diagram with ray incidence angles color coded?
	
DelayGraph = False      # produce an interactive spot diagram with delays in 3rd dimension?
IntensityGraph = True   # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
IncidenceGraph = False  # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?
	
DelayMirrorProjection = True      # produce a plot of the ray delays at the detector projected onto the mirror surface?
IntensityMirrorProjection = False # produce a plot of the ray intensities at the detector projected onto the mirror surface?
IncidenceMirrorProjection = False # produce a plot of the ray incidence angles at the detector projected onto the mirror surface?


#%%
""" a toroidal in 2f-2f config, but then additionally a normal-incidence spherical mirror to microfocus """

#% source properties
Divergence = 30e-3/2  # half-angle in rad
SourceSize = 0   # in mm
Wavelength = 50e-6  # 50 nm, only plays a role for the Airy spot size show in some plots as a reference
DeltaFT = 0.5       # 500 as, only appears as a reference in some plots
NumberRays = 1000

#% Mirror

Focal1 = 300
Focal2 = 50
AngleIncidence = 80 #in deg
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal1, AngleIncidence)
ToroidalMirror = mmirror.MirrorToroidal( OptimalMajorRadius, OptimalMinorRadius,msupp.SupportRectangle(120,30))
SphericalMirror = mmirror.MirrorSpherical( 2*Focal2,  msupp.SupportRound(12))
#SphericalMirror = mmirror.MirrorPlane( msupp.SupportRound(12))

#% creating the optical chain
OpticsList = [ToroidalMirror, SphericalMirror]
IncidenceAngleList = [AngleIncidence,0] #in deg
IncidencePlaneAngleList = [0,0]

#% detector
ReflectionNumber = -1 # analyse ray bundle after last optical element

ManualDetector = False # do not place the detector manually, i.e. let ART do it automatically
DistanceDetector = Focal2 # first set the detector at this distance from the last optical element
AutoDetectorDistance = True  # but then search for the optimal detector distance and shift the detector there


#%% LOOP OVER THE DISTANCE BETWEEN THE TWO MIRRORS
# initialize lists to contains all the results we want to keep from the loop iterations
DistanceS = [] 
RayListHistoryS = []
OpticalChainS = []
DetectorS = []
OptDistanceS = []
OptSizeSpotS = []
OptDurationS = []
DistanceBestSpotS = []
BestSizeS = []
DistanceBestDurationS = []
BestDurationS = []

#for Distance in tqdm(list(range(100, 500, 20))+list(range(600, 2000, 100))):
#for Distance in tqdm(list(range(50, 250, 50))+list(range(275, 2*Focal1-100, 25))):
#for Distance in tqdm(list(range(55, 75, 1))):
for Distance in [350]:
       
    DistanceList = [2*Focal1, Distance]
    
    OpticalChain = mp.OEPlacement(Divergence, SourceSize, NumberRays, OpticsList, DistanceList, IncidenceAngleList, IncidencePlaneAngleList)

    #firstMirror = OpticalChain[1]
    #firstMirror.position = firstMirror.position + 4*firstMirror.majoraxis
    
    #secondMirror = OpticalChain[-1]
    #secondMirror.position = secondMirror.position + 3*secondMirror.majoraxis
    
    #%%#####################################################################################################################
    """ Do THE ACTUAL RAYTRACING calculations """
    """ You can loop over this block to simulate a range of optical setups, e.g. a varying distance """
    
    exec( open("ART/DoART.py").read() )
    
    # append the lists to save the results of all iterations
    DistanceS.append(Distance)
    RayListHistoryS.append(RayListHistory)
    OpticalChainS.append(OpticalChain)
    DetectorS.append(Detector)
    OptDistanceS.append(OptDistance)
    OptSizeSpotS.append(OptSizeSpot)
    OptDurationS.append(OptDuration)
    DistanceBestSpotS.append(DistanceBestSpot)
    BestSizeS.append(BestSize)
    DistanceBestDurationS.append(DistanceBestDuration)
    BestDurationS.append(BestDuration)
    ########################################################################################################################

    # mlab.options.offscreen = True
    # rayRenderFig = mplots.RayRenderGraph(RayListHistory,OpticalChain,Detector.get_distance()*1.5 ,maxRaysToRender)
    # name = "f-{foo:04d}mm-f_with_f="+str(int(Focal))+"mm.png"
    # mlab.savefig(name.format(foo=int(Distance)), figure=rayRenderFig, magnification=3)

#%% create a pandas DataFrame to contain all the saved lists with descriptive column titles, which can then be saved to an efficient hdf file
resultsDF = pd.DataFrame({'Mirror Distance [mm]':DistanceS, 'RayListHistory':RayListHistoryS, 'OpticalChain':OpticalChainS, 'OptDetector':DetectorS, 'OptDetDistance [mm]':OptDistanceS, 'OptSpotSize [micron]':[x*1000 for x in OptSizeSpotS], 'OptDuration [fs]':OptDurationS, 'DistanceBestSpot [mm]':DistanceBestSpotS, 'BestSpotSize [micron]':[x*1000 for x in BestSizeS], 'DistanceBestDuration [mm]':DistanceBestDurationS, 'BestDuration [fs]':BestDurationS})
resultsDF.plot(kind='scatter',x='Mirror Distance [mm]',y='OptDuration [fs]')

#%%
# """possibly store results"""
# name = "f-x-f_for_"+str(int(Divergence*1000))+"mrad.h5"
# resultsDF.to_hdf(name, 'resultsDF'+str(int(Divergence*1000))+'mrad')


#%%#####################################################################################################################  
# Select the setup that gives the shortest duration at its best focus    
if len(resultsDF)>1:
    MinDurationRow = resultsDF['OptDuration [fs]'].idxmin()
    
    OpticalChain = resultsDF.at[MinDurationRow,'OpticalChain']
    RayListHistory = resultsDF.at[MinDurationRow,'RayListHistory']
    RayListAnalysed = RayListHistory[ReflectionNumber]
    Detector = resultsDF.at[MinDurationRow,'OptDetector']

#%
"""  Go through the BUILT-IN PLOT OPTIONS """

exec( open("ART/DoPLOTS.py").read() )
########################################################################################################################



