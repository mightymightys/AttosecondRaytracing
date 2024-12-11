"""
Module containing some useful functions for generating various plots.
It exists mostly to avoid cluttering the main plotting module.

Created in November 2024

@author: André Kalouguine + Stefan Haessler + Anthony Guillaume
"""
# %% Module imports
import weakref
import numpy as np
import colorcet as cc

import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleProcessing as mp


# %% Definition of observer class to make plots interactive
class Observable():
    """
    Observer class to make several plots interactive simultaneously.
    It encapsulates a single value.
    When it's modified, it notifies all registered observers.
    """
    def __init__(self, value):
        self._value = value
        self._observers = set()
        self._calculation = set()
    
    def register_calculation(self, callback):
        """
        Register a calculation to be performed when the value is modified before notifying the observers.
        For instance, this can be used to perform the ray tracing calculation when the value is modified.
        The other callbacks will be notified after the calculation is performed and will update the plots.
        """
        self._calculation.add(callback)

    def unregister_calculation(self, callback):
        self._calculation.discard(callback)
    
    @property
    def value(self):
        return self._value
    
    @value.setter
    def value(self, new_value):
        self._value = new_value
        for callback in self._calculation:
            callback(new_value)
        self.notify(new_value)
    
    def register(self, callback):
        self._observers.add(callback)

    def unregister(self, callback):
        self._observers.discard(callback)

    def notify(self, event):
        for callback in self._observers:
            callback(event)


# %% Utility functions
def generate_distinct_colors(num_colors):
    """
    Utility function to generate a list of distinct colors for plotting.

    Parameters
    ----------
        num_colors : int
            The number of colors to generate.

    Returns
    -------
        distinct_colors : list
            List of distinct colors.
    """
    # Get a color palette from colorcet
    palette = cc.glasbey

    # Make sure the number of colors does not exceed the palette length
    num_colors = min(num_colors, len(palette))

    # Slice the palette to get the desired number of colors
    distinct_colors = palette[:num_colors]

    return distinct_colors


def _getDetectorPoints(RayListAnalysed, Detector) -> tuple[np.ndarray, np.ndarray, float, float]:
    """
    Prepare the ray impact points on the detector in a format used for the plotting,
    and along the way also calculate the "spotsize" of this point-cloud on the detector.

    Parameters
    ----------
        RayListAnalysed : list(Ray)
            A list of objects of the ModuleOpticalRay.Ray-class.

        Detector : Detector
            An object of the ModuleDetector.Detector-class.

    Returns
    -------
        DectectorPoint2D_Xcoord : np.ndarray with x-coordinates

        DectectorPoint2D_Ycoord : np.ndarray with y-coordinates

        FocalSpotSize : float

        SpotSizeSD : float
    """

    Points2D = Detector.get_2D_points(RayListAnalysed)
    Points2D -= np.mean(Points2D, axis=1)  # Centering the points
    X = Points2D[0][:,0] * 1e3 # To convert to µm
    Y = Points2D[0][:,1] * 1e3 # To convert to µm

    FocalSpotSize = float(mgeo.DiameterPointArray(Points2D[0]))
    SpotSizeSD = mp.StandardDeviation(Points2D[0])
    return X, Y, FocalSpotSize, SpotSizeSD

def show():
    plt.show(block=False)