# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 13:54:51 2023

@author: Haessler
"""
DefaultAnalysisOptions = {
    "verbose": True,  # print intermediate results and info in the console?
    "plot_Render": False,  # render optical elements and rays, and how many rays to render?
    "maxRaysToRender": 200,
    "OEPointsToRender": 3000,
    "OEPointsScale": 5,
    "draw_mesh": False,
    "cycle_ray_colors": False,
    "DrawAiryAndFourier": True,  # Draw Airy spot and Fourier-limited duration in the following plots?
    "plot_SpotDiagram": False,  # produce an interactive spot diagram without color coding the spots?
    "plot_DelaySpotDiagram": False,  # produce an interactive spot diagram with ray delays color coded?
    "plot_IntensitySpotDiagram": False,  # produce an interactive spot diagram with ray intensities color coded?
    "plot_IncidenceSpotDiagram": False,  # produce an interactive spot diagram with ray incidence angles color coded?
    "plot_DelayGraph": False,  # produce an interactive spot diagram with delays in 3rd dimension?
    "plot_IntensityGraph": False,  # produce an interactive spot diagram with delays in 3rd dimension and ray intensities color coded?
    "plot_IncidenceGraph": False,  # produce an interactive spot diagram with delays in 3rd dimension and ray incidence angles color coded?
    "plot_DelayMirrorProjection": False,  # produce a plot of the ray delays at the detector projected onto the mirror surface?
    "plot_IntensityMirrorProjection": False,  # produce a plot of the ray intensities at the detector projected onto the mirror surface?
    "plot_IncidenceMirrorProjection": False,  # produce a plot of the ray incidence angles at the detector projected onto the mirror surface?
    "save_results": True,  # save the simulation results to disk, to analyse later
}

# %%
DefaultSourceProperties = {
    "Divergence": 0,  # half-angle in rad
    "SourceSize": 0,  # diameter in mm
    "Wavelength": 50e-6,  # 50 nm, typical XUV value; only plays a role for the Airy spot size show in some plots as a reference
    "DeltaFT": 1,  # 1 fs, because why not
    "NumberRays": 1000,  # launch 1000 rays in the beginning
}

# %%
DefaultDetectorOptions = {
    "ReflectionNumber": -1,  # analyse ray bundle after last optical element
    "ManualDetector": False,  # do not place the detector manually, i.e. let ART do it automatically
    "DetectorCentre": None,  # in case of manual detector
    "DetectorNormal": None,  # in case of manual detector
    "DistanceDetector": None,  # set the detector at this distance from the last optical element
    "AutoDetectorDistance": False,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    "OptFor": "intensity",  # metric for which to optimize the detector position
}
