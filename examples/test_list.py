# -*- coding: utf-8 -*-
"""
Created in Apr 2020

@author: Stefan Haessler
"""
#%% Modules
#import copy

import numpy as np
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleMask as mmask
import ARTcore.ModuleSource as mos
import ARTcore.ModuleOpticalChain as moc
import ART.ModuleAnalysisAndPlots as maap
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleDetector as mdet
from ART.ARTmain import run_ART
from copy import copy
import matplotlib.pyplot as plt
from scipy.stats import linregress
import ART.ModuleAnalysis as man
import time
start_time = time.time()


#%%########################################################################
Spectrum = mos.UniformSpectrum(lambdaMin=30e-6, lambdaMax=800e-6)
#Spectrum = mos.SingleWavelengthSpectrum(800e-6)
PowerDistribution = mos.GaussianPowerDistribution(1, 2, 50e-3)
Positions = mos.PointRayOriginsDistribution(mgeo.Origin)
Directions = mos.ConeRayDirectionsDistribution(mgeo.Vector([1,0,0]), 50e-3)
Source = mos.SimpleSource(Spectrum, PowerDistribution, Positions, Directions)

ChainDescription = "2 toroidal mirrors in f-d-f config, i.e. approx. collimation, propagation, and the refocus "



# %% Define the optical elements
SupportMask = msupp.SupportRoundHole(Radius=20, RadiusHole=14/2, CenterHoleX=0, CenterHoleY=0) 
Mask = mmask.Mask(SupportMask)
MaskSettings = {
    'OpticalElement' : Mask,
    'Distance' : 400,
    'IncidenceAngle' : 0,
    'IncidencePlaneAngle' : 0,
    'Description' : "Mask for selecting rays",
    'Alignment' : 'support_normal',
}

Focal = 500
AngleIncidence = 80 #in deg
OptimalMajorRadius, OptimalMinorRadius = mmirror.ReturnOptimalToroidalRadii(Focal, AngleIncidence)
SupportToroidal = msupp.SupportRectangle(150, 32)

ToroidalMirrorA = mmirror.MirrorToroidal(SupportToroidal, OptimalMajorRadius, OptimalMinorRadius)
ToroidalASettings = {
    'OpticalElement' : ToroidalMirrorA,
    'Distance' : Focal-MaskSettings['Distance'],
    'IncidenceAngle' : 0,
    'IncidencePlaneAngle' : 0,
    'Description' : "First parabola for collimation",
}

ToroidalMirrorB = mmirror.MirrorToroidal(SupportToroidal,OptimalMajorRadius, OptimalMinorRadius)
ToroidalBSettings = {
    'OpticalElement' : ToroidalMirrorB,
    'Distance' : None,
    'IncidenceAngle' : 0,
    'IncidencePlaneAngle' : 0,
    'Description' : "First parabola for collimation",
}

Det = mdet.InfiniteDetector()
Detectors = {
    "Focus": (Det, -1) # -1 means that the detector is placed at the last optical element
}


Distances = np.linspace(Focal-200, Focal+200, 100)
FocalSizes = []

for d in Distances:
    toroidalBSettings = copy(ToroidalBSettings)
    toroidalBSettings['OpticalElement'] = copy(ToroidalMirrorB)
    toroidalBSettings['Distance'] = d
    toroidalASettings = copy(ToroidalBSettings)
    toroidalASettings['OpticalElement'] = copy(ToroidalMirrorB)
    maskSettings = copy(MaskSettings)
    maskSettings['OpticalElement'] = copy(Mask)
    OpticsList = [maskSettings,toroidalASettings, toroidalBSettings]
    AlignedOpticalElements = mp.OEPlacement(OpticsList)
    print(ToroidalBSettings["OpticalElement"].basis)
    AlignedOpticalChain = moc.OpticalChain(Source(1000), AlignedOpticalElements, Detectors, ChainDescription)
    rays= AlignedOpticalChain.get_output_rays()
    Det.autoplace(rays[-1], 390)
    Det.optimise_distance(AlignedOpticalChain.get_output_rays()[-1], [200,600], Det._spot_size, maxiter=10, tol=1e-14)
    FocalSizes.append(Det.distance)


# AlignedOpticalChain.rotate_OE(1, "localnormal", "roll", 180)
# AlignedOpticalChain.partial_realign(2,3, DistanceList[2:3], IncidenceAngleList[2:3], IncidencePlaneAngleList[2:3])

#Detector = setup_detector(AlignedOpticalChain, DetectorOptions, AlignedOpticalChain.get_output_rays()[-1])
#maap.SpotDiagram(AlignedOpticalChain.get_output_rays()[-1], Detector)
#maap.SpotDiagram(AlignedOpticalChain, ColorCoded="Incidence", DrawAiryAndFourier=True)
#maap.RayRenderGraph(AlignedOpticalChain, EndDistance=500, OEpoints=5000, cycle_ray_colors=True, impact_points=True, DetectedRays=True)
#X,Y,Z = man.get_planewavefocus(AlignedOpticalChain, DetectorName="Focus", size=None, Nrays=1000, resolution=100)
#plt.pcolormesh(X*1e3,Y*1e3,Z)
#plt.show()