# -*- coding: utf-8 -*-
"""
Created in Jan 2023

@author: Haessler & Semptum
"""
# %% Modules
# import copy
import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp
import ART.ModuleDefects as mdef
import numpy as np
from ARTmain import main

# %%########################################################################
""" Source properties """

SourceProperties = {
    "Divergence": 0,  # half-angle in rad, 0 for a plane wave
    "SourceSize": 30,  # diameter in mm, 0 for a point source
    "Wavelength": 800e-6,  # 800 nm
    "DeltaFT": 0,  # in fs
    "NumberRays": 1000,  # launch 1000 rays in the beginning
}


# %%
""" OPTICAL SETUP """
Description = "TODO"

# Mirror
Support = msupp.SupportRound(20)
# Support =  msupp.SupportRound(50)

FocalEffective = 25.4  # in mm

Mirror = mmirror.MirrorParabolic(FocalEffective, 0, Support)

RMS_values = np.exp(np.linspace(-9,0, num=100)*np.log(10))
Defects = [mdef.Fourrier(Support, RMS = i, smallest=1) for i in RMS_values]
DeformedMirrors = [mmirror.DeformedMirror(Mirror, [i]) for i in Defects]
STDs = dict()
DeformedMirrors = [DeformedMirrors[0]]

for DeformedMirror in DeformedMirrors:
    OpticsList = [DeformedMirror]  # , DeformedEllipse2]
    DistanceList = [15]  # 50.8,
    IncidenceAngleList = [10]
    
    
    OpticalChain = mp.OEPlacement(SourceProperties, OpticsList, DistanceList, IncidenceAngleList, Description=Description)

    DetectorOptions = {
        "ReflectionNumber": -1,
        "ManualDetector": False,
        "DistanceDetector": FocalEffective,
        "AutoDetectorDistance": False,
        "OptFor": "intensity",
    }
    
    AnalysisOptions = {
        "verbose": True,
        "plot_Render": True,
        "maxRaysToRender": 150,
        "DrawAiryAndFourier": False,
        "plot_SpotDiagram": False,
        "plot_DelaySpotDiagram": False,
        "plot_IntensitySpotDiagram": False,
        "plot_IncidenceSpotDiagram": False,
        "plot_DelayGraph": False,
        "plot_IntensityGraph": False,
        "plot_IncidenceGraph": False,
        "plot_DelayMirrorProjection": False,
        "plot_IntensityMirrorProjection": False,
        "plot_IncidenceMirrorProjection": False,
        "save_results": False,
        "OEPointsToRender": 4000,
        "OEPointsScale": 0,
    }

    kept_data = main(OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions)
    STDs[DeformedMirror.DeformationList[0].RMS()] = kept_data["SpotSizeSD"]

RMS_vals = sorted(list(STDs.keys()))
stds = [STDs[i] for i in RMS_vals]