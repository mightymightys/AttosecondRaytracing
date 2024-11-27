"""
Provides functions for analysing the alignment tolerance of a system.

Mostly it consists in mis-aligning the optical elements of the system and checking the impact on the focal spot and delays.

Created in Nov 2024

@author: Andre Kalouguine
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


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
