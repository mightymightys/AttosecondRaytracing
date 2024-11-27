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
from ART.ARTmain import run_ART
from copy import copy
import matplotlib.pyplot as plt
from scipy.stats import linregress
from ART.ARTmain import setup_detector

#%%########################################################################
""" Source properties """

SourceProperties = {
    'Divergence' : 80e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 50e-6, # 40 nm
    'DeltaFT'    : 0.3,   # in fs
    'NumberRays' : 5000   # launch 5000 rays in the beginning 
}

DetectorOptions = {
    'ReflectionNumber' : -1, # analyse ray bundle after last optical element
    'ManualDetector' : False, # do not place the detector manually, i.e. let ART do it automatically
    'DistanceDetector' : 0 ,  # set the detector at this distance from the last optical element
    'AutoDetectorDistance' : False,  # search for the optimal detector distance and shift the detector there to use it for further analysis ?
    'OptFor' : "intensity"   # metric for which to optimize the detector position
}

#%%########################################################################

SourceRays = mos.PointSource(
    np.zeros(3),
    np.array([1,0,0]),
    SourceProperties["Divergence"],
    SourceProperties["NumberRays"], 
    SourceProperties["Wavelength"])

Source = mos.ApplyGaussianIntensityToRayList(
        SourceRays, 1 / np.e**2
    )

Description = "2 equal large-off-axis-angle parabolas for collimation and refocusing "

Support = msupp.SupportRectangle(40,40)
offAxisAngle = 160 #in deg
FocalEffective = 600 # in mm
Parabola = mmirror.MirrorParabolic(FocalEffective, offAxisAngle, Support) 
ParabolaBis = mmirror.MirrorParabolic(FocalEffective, offAxisAngle, Support) 

SupportPlane = msupp.SupportRound(75)
PlaneMirror = mmirror.MirrorPlane(SupportPlane) 
PlaneMirrorB = mmirror.MirrorPlane(SupportPlane) 

maskdistance = 390

RadiusHole = 50e-3/2*maskdistance
RadiusMask = 50
SupportMask = msupp.SupportRoundHole(Radius=RadiusMask, RadiusHole = RadiusHole, CenterHoleX=0, CenterHoleY=0) 
Mask = mmask.Mask(SupportMask)



OpticsList = [Mask, Parabola, PlaneMirror, ParabolaBis]
DistanceList = [maskdistance, FocalEffective-maskdistance, 300, 300]
IncidenceAngleList = [0, offAxisAngle, 70, 0]
IncidencePlaneAngleList = [0, 0, 180, 0]

AlignedOpticalElements = mp.OEPlacement(OpticsList, DistanceList,IncidenceAngleList,IncidencePlaneAngleList)
AlignedOpticalChain = moc.OpticalChain(Source, AlignedOpticalElements, Description)
# AlignedOpticalChain.rotate_OE(1, "localnormal", "roll", 180)
# AlignedOpticalChain.partial_realign(2,3, DistanceList[2:3], IncidenceAngleList[2:3], IncidencePlaneAngleList[2:3])

Detector = setup_detector(AlignedOpticalChain, DetectorOptions, AlignedOpticalChain.get_output_rays()[0])
maap.SpotDiagram(AlignedOpticalChain.get_output_rays()[-1], Detector)

def fit_and_plot(ax, displacements, results, movement_label, index):
    # Convert data to numpy arrays for handling
    displacements = np.array(displacements)
    results = np.array(results)

    # Handling the symmetry: use absolute displacements for fitting
    fit_results = linregress(np.abs(displacements), results)
    slope, intercept, r_value, p_value, std_err = fit_results
    
    # Prepare data for plotting the fit line
    fit_line_x = np.linspace(np.min(np.abs(displacements)), np.max(np.abs(displacements)), 100)
    fit_line_y = intercept + slope * fit_line_x
    
    # Scatter plot of the original data
    ax.scatter(displacements, results, label=f'{movement_label.capitalize()} (Data)')
    
    # Plot the fit line
    ax.plot(fit_line_x, fit_line_y, 'r-', label=f'{movement_label.capitalize()} Fit (slope={slope:.2f}±{std_err:.2f})')
    
    ax.legend()
    ax.set_xlabel("Movement [units]")
    ax.set_ylabel("Pulse duration SD [fs]")
    ax.set_title(f"Effect of {movement_label} the optical element at index {index}")


def vary_orientation(optical_chain, index, reference="out", num_simulations=100, rotation_range=0.1):
    """ Varies the orientation of a specified optical element and plots the effect on pulse duration."""
    results = {'pitch': [], 'roll': [], 'yaw': []}
    angles = {'pitch': [], 'roll': [], 'yaw': []}
    
    for rotation in ['pitch', 'roll', 'yaw']:
        for _ in range(num_simulations):
            misaligned_optical_elements = copy(optical_chain.optical_elements)
            misaligned_optical_chain = moc.OpticalChain(optical_chain.source_rays, misaligned_optical_elements, "Misaligned "+Description)
            angle = (np.random.rand() - 0.5) * rotation_range * 2
            misaligned_optical_chain.rotate_OE(index, reference, rotation, angle)
            SD = run_ART(misaligned_optical_chain, SourceProperties, DetectorOptions, AnalysisOptions)[-1]
            results[rotation].append(SD)
            angles[rotation].append(angle)
    fig, ax =  plt.subplots()
    for rotation in ['pitch', 'roll', 'yaw']:
        fit_and_plot(ax, angles[rotation], results[rotation], rotation, index)
    
    ax.legend()
    ax.set_xlabel("Rotation [°]")
    ax.set_ylabel("Pulse duration SD [fs]")
    ax.set_title(f"Effect of rotating the optical element at index {index}")
    plt.show()
    return fig

def vary_position(optical_chain, index, reference="out", num_simulations=100, displacement_range=1.0):
    """ Varies the position of a specified optical element and plots the effect on pulse duration."""
    results = {'along': [], 'in_plane': [], 'out_plane': []}
    displacements = {'along': [], 'in_plane': [], 'out_plane': []}
    

    for direction in ['along', 'in_plane', 'out_plane']:
        for _ in range(num_simulations):
            misaligned_optical_elements = copy(optical_chain.optical_elements)
            misaligned_optical_chain = moc.OpticalChain(optical_chain.source_rays, misaligned_optical_elements, "Misaligned "+Description)
            displacement = (np.random.rand() - 0.5) * displacement_range * 2 
            misaligned_optical_chain.shift_OE(index, reference, direction, displacement)
            SD = run_ART(misaligned_optical_chain, SourceProperties, DetectorOptions, AnalysisOptions)[-1]
            results[direction].append(SD)
            displacements[direction].append(displacement)
    fig, ax =  plt.subplots()
    for direction in ['along', 'in_plane', 'out_plane']:
        fit_and_plot(ax, displacements[direction], results[direction], direction.replace('_', ' ').capitalize(), index)
    
    ax.legend()
    ax.set_xlabel("Displacement [mm]")
    ax.set_ylabel("Pulse duration SD [fs]")
    ax.set_title(f"Effect of translating the optical element at index {index}")
    plt.show()
    return fig
