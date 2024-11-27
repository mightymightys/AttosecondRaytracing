"""
Provides functions for analysis of the output ray-bundles calculated through ray-tracing, and also for the ART's standard visualization options.

Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler + Andre Kalouguine
"""
# %% Module imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
import pyvistaqt as pvqt
import colorcet as cc
from colorsys import rgb_to_hls, hls_to_rgb
import logging

logger = logging.getLogger(__name__)


import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo
import ARTcore.ModuleDetector as mdet
import ART.ModulePlottingMethods as mpm
import ART.ModulePlottingUtilities as mpu
from ART.ModulePlottingUtilities import Observable

import ARTcore.ModuleOpticalChain as moc
import ART.ModuleAnalysis as man
import itertools
from copy import copy


# %% Spot diagram on detector
def SpotDiagram(OpticalChain, DetectorName = "Focus", 
                DrawAiryAndFourier=False, 
                DrawFocalContour=False,
                DrawFocal=False,
                ColorCoded=None,
                Observer = None) -> plt.Figure:
    """
    Produce an interactive figure with the spot diagram on the selected Detector.
    The detector distance can be shifted with the left-right cursor keys. Doing so will actually move the detector.
    If DrawAiryAndFourier is True, a circle with the Airy-spot-size will be shown.
    If DrawFocalContour is True, the focal contour calculated from some of the rays will be shown.
    If DrawFocal is True, a heatmap calculated from some of the rays will be shown. 
    The 'spots' can optionally be color-coded by specifying ColorCoded, which can be one of ["Intensity","Incidence","Delay"].

    Parameters
    ----------
        RayListAnalysed : list(Ray)
            List of objects of the ModuleOpticalRay.Ray-class.

        Detector : Detector
            An object of the ModuleDetector.Detector-class.

        DrawAiryAndFourier : bool, optional
            Whether to draw a circle with the Airy-spot-size. The default is False.
        
        DrawFocalContour : bool, optional
            Whether to draw the focal contour. The default is False.

        DrawFocal : bool, optional
            Whether to draw the focal heatmap. The default is False.

        ColorCoded : str, optional
            Color-code the spots according to one of ["Intensity","Incidence","Delay"]. The default is None.

        Observer : Observer, optional
            An observer object. If none, then we just create a copy of the detector and move it when pressing left-right. 
            However, if an observer is specified, then we will change the value of the observer and it will issue 
            the required callbacks to update several plots at the same time.

    Returns
    -------
        fig : matlplotlib-figure-handle.
            Shows the interactive figure.
    """
    Detector, Index = OpticalChain.detectors[DetectorName]
    movingDetector = copy(Detector) # We will move this detector when pressing left-right
    if Observer is None:
        detectorPosition = Observable(movingDetector.distance) # We will observe the distance of this detector
    else:
        detectorPosition = Observer
        movingDetector.distance = detectorPosition.value

    detectorPosition.register_calculation(lambda x: movingDetector.set_distance(x))

    RayListAnalysed = OpticalChain.get_output_rays()[Index]

    NumericalAperture = man.ReturnNumericalAperture(RayListAnalysed, 1)  # NA determined from final ray bundle
    MaxWavelength = np.max([i.wavelength for i in RayListAnalysed])
    if DrawAiryAndFourier:
        AiryRadius = man.ReturnAiryRadius(MaxWavelength, NumericalAperture) * 1e3  # in µm
    else:
        AiryRadius = 0
    
    if DrawFocalContour or DrawFocal:
        X,Y,Z = man.get_diffractionfocus(OpticalChain, movingDetector, Index)
        Z/=np.max(Z) 

    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = mpu._getDetectorPoints(
        RayListAnalysed, movingDetector
    )

    match ColorCoded:
        case "Intensity":
            IntensityList = [k.intensity for k in RayListAnalysed]
            z = np.asarray(IntensityList)
            zlabel = "Intensity (arb.u.)"
            title = "Intensity + Spot Diagram\n press left/right to move detector position"
            addLine = ""
        case "Incidence":
            IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed]  # degree
            z = np.asarray(IncidenceList)
            zlabel = "Incidence angle (deg)"
            title = "Ray Incidence + Spot Diagram\n press left/right to move detector position"
            addLine = ""
        case "Delay":
            DelayList = movingDetector.get_Delays(RayListAnalysed)
            DurationSD = mp.StandardDeviation(DelayList)
            z = np.asarray(DelayList)
            zlabel = "Delay (fs)"
            title = "Delay + Spot Diagram\n press left/right to move detector position"
            addLine = "\n" + "{:.2f}".format(DurationSD) + " fs SD"
        case _:
            z = "red"
            title = "Spot Diagram\n press left/right to move detector position"
            addLine = ""

    distStep = min(50, max(0.0005, round(FocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000))  # in mm

    plt.ion()
    fig, ax = plt.subplots()
    if DrawFocal:
        focal = ax.pcolormesh(X*1e3,Y*1e3,Z)
    if DrawFocalContour:
        levels = [1/np.e**2, 0.5]
        contour = ax.contourf(X*1e3, Y*1e3, Z, levels=levels, cmap='gray')

    if DrawAiryAndFourier:
        theta = np.linspace(0, 2 * np.pi, 100)
        x = AiryRadius * np.cos(theta)
        y = AiryRadius * np.sin(theta)  #
        ax.plot(x, y, c="black")
        

    foo = ax.scatter(
        DectectorPoint2D_Xcoord,
        DectectorPoint2D_Ycoord,
        c=z,
        s=15,
        label="{:.3f}".format(detectorPosition.value) + " mm\n" + "{:.1f}".format(SpotSizeSD * 1e3) + " \u03BCm SD" + addLine,
    )

    axisLim = 1.1 * max(AiryRadius, 0.5 * FocalSpotSize * 1000)
    ax.set_xlim(-axisLim, axisLim)
    ax.set_ylim(-axisLim, axisLim)

    if ColorCoded == "Intensity" or ColorCoded == "Incidence" or ColorCoded == "Delay":
        cbar = fig.colorbar(foo)
        cbar.set_label(zlabel)

    ax.legend(loc="upper right")
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_title(title)
    # ax.margins(x=0)


    def update_plot(new_value):
        nonlocal movingDetector, ColorCoded, zlabel, cbar, detectorPosition, foo, distStep, focal, contour, levels, Index, RayListAnalysed

        newDectectorPoint2D_Xcoord, newDectectorPoint2D_Ycoord, newFocalSpotSize, newSpotSizeSD = mpu._getDetectorPoints(
            RayListAnalysed, movingDetector
        )

        if DrawFocal:
            focal.set_array(Z)
        if DrawFocalContour:
            levels = [1/np.e**2, 0.5]
            for coll in contour.collections:
                coll.remove()  # Remove old contour lines
            contour = ax.contourf(X * 1e3, Y * 1e3, Z, levels=levels, cmap='gray')
        
        xy = foo.get_offsets()
        xy[:, 0] = newDectectorPoint2D_Xcoord
        xy[:, 1] = newDectectorPoint2D_Ycoord
        foo.set_offsets(xy)


        if ColorCoded == "Delay":
            newDelayList = np.asarray(movingDetector.get_Delays(RayListAnalysed))
            newDurationSD = mp.StandardDeviation(newDelayList)
            newaddLine = "\n" + "{:.2f}".format(newDurationSD) + " fs SD"
            foo.set_array(newDelayList)
            foo.set_clim(min(newDelayList), max(newDelayList))
            cbar.update_normal(foo)
        else:
            newaddLine = ""

        foo.set_label(
            "{:.3f}".format(detectorPosition.value) + " mm\n" + "{:.1f}".format(newSpotSizeSD * 1e3) + " \u03BCm SD" + newaddLine
        )
        ax.legend(loc="upper right")

        axisLim = 1.1 * max(AiryRadius, 0.5 * newFocalSpotSize * 1000)
        ax.set_xlim(-axisLim, axisLim)
        ax.set_ylim(-axisLim, axisLim)

        distStep = min(
            50, max(0.0005, round(newFocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000)
        )  # in mm

        fig.canvas.draw_idle()


    def press(event):
        nonlocal detectorPosition, distStep
        if event.key == "right":
            detectorPosition.value += distStep
        elif event.key == "left":
            if detectorPosition.value > 1.5 * distStep:
                detectorPosition.value -= distStep
            else:
                detectorPosition.value = 0.5 * distStep
        else:
            return None

    fig.canvas.mpl_connect("key_press_event", press)

    plt.show()

    detectorPosition.register(update_plot)


    return fig, detectorPosition


# %% Delay graph on detector (3D spot diagram)
def _drawDelayGraph(RayListAnalysed, Detector, Distance, DeltaFT, DrawAiryAndFourier=False, ColorCoded=None, fig=None):
    """
    Draws the 3D-delay-spot-diagram for a fixed detector. See more doc in the function below.
    """
    NumericalAperture = man.ReturnNumericalAperture(RayListAnalysed, 1)  # NA determined from final ray bundle
    Wavelength = RayListAnalysed[0].wavelength
    AiryRadius = man.ReturnAiryRadius(Wavelength, NumericalAperture) * 1e3  # in µm

    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = mpu._getDetectorPoints(
        RayListAnalysed, Detector
    )

    DelayList = Detector.get_Delays(RayListAnalysed)
    DurationSD = mp.StandardDeviation(DelayList)

    if ColorCoded == "Intensity":
        IntensityList = [k.intensity for k in RayListAnalysed]
    elif ColorCoded == "Incidence":
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed]  # in degree

    plt.ion()
    if fig is None:
        fig = plt.figure()
    else:
        fig.clear(keep_observers=True)

    ax = Axes3D(fig)
    fig.add_axes(ax)
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")
    ax.set_zlabel("Delay (fs)")

    labelLine = (
        "{:.3f}".format(Distance)
        + " mm\n"
        + "{:.1f}".format(SpotSizeSD * 1e3)
        + " \u03BCm SD\n"
        + "{:.2f}".format(DurationSD)
        + " fs SD"
    )
    if ColorCoded == "Intensity":
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, DelayList, s=4, c=IntensityList, label=labelLine)
        # plt.title('Delay + Intensity graph\n press left/right to move detector position')
        ax.set_title("Delay + Intensity graph\n press left/right to move detector position")

    elif ColorCoded == "Incidence":
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, DelayList, s=4, c=IncidenceList, label=labelLine)
        # plt.title('Delay + Incidence graph\n press left/right to move detector position')
        ax.set_title("Delay + Incidence graph\n press left/right to move detector position")
    else:
        ax.scatter(DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, DelayList, s=4, c=DelayList, label=labelLine)
        # plt.title('Delay graph\n press left/right to move detector position')
        ax.set_title("Delay graph\n press left/right to move detector position")

    ax.legend(loc="upper right")

    if DrawAiryAndFourier:
        x = np.linspace(-AiryRadius, AiryRadius, 40)
        z = np.linspace(np.mean(DelayList) - DeltaFT * 0.5, np.mean(DelayList) + DeltaFT * 0.5, 40)
        x, z = np.meshgrid(x, z)
        y = np.sqrt(AiryRadius**2 - x**2)
        ax.plot_wireframe(x, y, z, color="grey", alpha=0.1)
        ax.plot_wireframe(x, -y, z, color="grey", alpha=0.1)

    axisLim = 1.1 * max(AiryRadius, 0.5 * FocalSpotSize * 1000)
    ax.set_xlim(-axisLim, axisLim)
    ax.set_ylim(-axisLim, axisLim)
    # ax.set_zlim(-axisLim/3*10, axisLim/3*10) #same scaling as spatial axes 
    
    plt.show()
    fig.canvas.draw()

    return fig, NumericalAperture, AiryRadius, FocalSpotSize

def DelayGraph(OpticalChain, DetectorName, DeltaFT: (int, float), 
        DrawAiryAndFourier=False, 
        ColorCoded=None,
        Observer = None
    ) -> plt.Figure:
    """
    Produce a an interactive figure with a spot diagram resulting from the RayListAnalysed
    hitting the Detector, with the ray-delays shown in the 3rd dimension.
    The detector distance can be shifted with the left-right cursor keys.
    If DrawAiryAndFourier is True, a cylinder is shown whose diameter is the Airy-spot-size and
    whose height is the Fourier-limited pulse duration 'given by 'DeltaFT'.
    
    The 'spots' can optionally be color-coded by specifying ColorCoded as ["Intensity","Incidence"].

    Parameters
    ----------
        RayListAnalysed : list(Ray)
            List of objects of the ModuleOpticalRay.Ray-class.

        Detector : Detector
            An object of the ModuleDetector.Detector-class.

        DeltaFT : (int, float)
            The Fourier-limited pulse duration. Just used as a reference to compare the temporal spread
            induced by the ray-delays.

        DrawAiryAndFourier : bool, optional
            Whether to draw a cylinder showing the Airy-spot-size and Fourier-limited-duration.
            The default is False.

        ColorCoded : str, optional
            Color-code the spots according to one of ["Intensity","Incidence"].
            The default is None.

    Returns
    -------
        fig : matlplotlib-figure-handle.
            Shows the interactive figure.
    """
    Det, Index = OpticalChain.detectors[DetectorName]
    Detector = copy(Det)
    if Observer is None:
        detectorPosition = Observable(Detector.distance)
    else:
        detectorPosition = Observer
        Detector.distance = detectorPosition.value
    
    RayListAnalysed = OpticalChain.get_output_rays()[Index]
    fig, NumericalAperture, AiryRadius, FocalSpotSize = _drawDelayGraph(
        RayListAnalysed, Detector, detectorPosition.value, DeltaFT, DrawAiryAndFourier, ColorCoded
    )

    distStep = min(50, max(0.0005, round(FocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000))  # in mm

    movingDetector = copy(Detector)

    def update_plot(new_value):
        nonlocal movingDetector, ColorCoded, detectorPosition, distStep, fig
        ax = fig.axes[0]
        cam = [ax.azim, ax.elev, ax._dist]
        fig, sameNumericalAperture, sameAiryRadius, newFocalSpotSize = _drawDelayGraph(
            RayListAnalysed, movingDetector, detectorPosition.value, DeltaFT, DrawAiryAndFourier, ColorCoded, fig
        )
        ax = fig.axes[0]
        ax.azim, ax.elev, ax._dist = cam
        distStep = min(
            50, max(0.0005, round(newFocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000)
        )
        return fig

    def press(event):
        nonlocal detectorPosition, distStep, movingDetector, fig
        if event.key == "right":
            detectorPosition.value += distStep
        elif event.key == "left":
            if detectorPosition.value > 1.5 * distStep:
                detectorPosition.value -= distStep
            else:
                detectorPosition.value = 0.5 * distStep

    fig.canvas.mpl_connect("key_press_event", press)
    detectorPosition.register(update_plot)
    detectorPosition.register_calculation(lambda x: movingDetector.set_distance(x))

    return fig, Observable


# %% Mirror projection
def MirrorProjection(OpticalChain, ReflectionNumber: int, ColorCoded=None, DetectorName="") -> plt.Figure:
    """
    Produce a plot of the ray impact points on the optical element with index 'ReflectionNumber'.
    The points can be color-coded according ["Incidence","Intensity","Delay"], where the ray delay is
    measured at the Detector.

    Parameters
    ----------
        OpticalChain : OpticalChain
           List of objects of the ModuleOpticalOpticalChain.OpticalChain-class.

        ReflectionNumber : int
            Index specifying the optical element on which you want to see the impact points.

        Detector : Detector, optional
            Object of the ModuleDetector.Detector-class. Only necessary to project delays. The default is None.

        ColorCoded : str, optional
            Specifies which ray property to color-code: ["Incidence","Intensity","Delay"]. The default is None.

    Returns
    -------
        fig : matlplotlib-figure-handle.
            Shows the figure.
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    Position = OpticalChain[ReflectionNumber].position
    q = OpticalChain[ReflectionNumber].orientation
    # n = OpticalChain.optical_elements[ReflectionNumber].normal
    # m = OpticalChain.optical_elements[ReflectionNumber].majoraxis

    RayListAnalysed = OpticalChain.get_output_rays()[ReflectionNumber]
    # transform rays into the mirror-support reference frame
    # (same as mirror frame but without the shift by mirror-centre)
    r0 = OpticalChain[ReflectionNumber].r0
    RayList = [r.to_basis(*OpticalChain[ReflectionNumber].basis) for r in RayListAnalysed]

    x = np.asarray([k.point[0] for k in RayList]) - r0[0]
    y = np.asarray([k.point[1] for k in RayList]) - r0[1]
    if ColorCoded == "Intensity":
        IntensityList = [k.intensity for k in RayListAnalysed]
        z = np.asarray(IntensityList)
        zlabel = "Intensity (arb.u.)"
        title = "Ray intensity projected on mirror              "
    elif ColorCoded == "Incidence":
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed]  # in degree
        z = np.asarray(IncidenceList)
        zlabel = "Incidence angle (deg)"
        title = "Ray incidence projected on mirror              "
    elif ColorCoded == "Delay":
        if OpticalChain.Detector is not None:
            z = np.asarray(Detector.get_Delays(RayListAnalysed))
            zlabel = "Delay (fs)"
            title = "Ray delay at detector projected on mirror              "
        else:
            raise ValueError("If you want to project ray delays, you must specify a detector.")
    else:
        z = "red"
        title = "Ray impact points projected on mirror"

    plt.ion()
    fig = plt.figure()
    ax = OpticalChain.optical_elements[ReflectionNumber].support._ContourSupport(fig)
    p = plt.scatter(x, y, c=z, s=15)
    if ColorCoded == "Delay" or ColorCoded == "Incidence" or ColorCoded == "Intensity":
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(p, cax=cax)
        cbar.set_label(zlabel)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.title(title, loc="right")
    plt.tight_layout()

    bbox = ax.get_position()
    bbox.set_points(bbox.get_points() - np.array([[0.01, 0], [0.01, 0]]))
    ax.set_position(bbox)
    plt.show()

    return fig


# %% Setup rendering
def _RenderRays(RayListHistory, EndDistance, maxRays=150):
    """
    Generates a list of Pyvista-meshes representing the rays in RayListHistory.
    This can then be employed to render the rays in a Pyvista-figure.

    Parameters
    ----------
        RayListHistory : list(list(Ray))
            A list of lists of objects of the ModuleOpticalRay.Ray-class.

        EndDistance : float
            The rays of the last ray bundle are drawn with a length given by EndDistance (in mm).

        maxRays : int, optional
            The maximum number of rays to render. Rendering all the traced rays is a insufferable resource hog
            and not required for a nice image. Default is 150.
    
    Returns
    -------
        meshes : list(pvPolyData)
            List of Pyvista PolyData objects representing the rays
    """
    meshes = []
    # Ray display
    for k in range(len(RayListHistory)):
        x = []
        y = []
        z = []
        if k != len(RayListHistory) - 1:
            knums = list(
                map(lambda x: x.number, RayListHistory[k])
            )  # make a list of all ray numbers that are still in the game
            if len(RayListHistory[k + 1]) > maxRays:
                rays_to_render = np.random.choice(RayListHistory[k + 1], maxRays, replace=False)
            else:
                rays_to_render = RayListHistory[k + 1]

            for j in rays_to_render:
                indx = knums.index(j.number)
                i = RayListHistory[k][indx]
                Point1 = i.point
                Point2 = j.point
                x += [Point1[0], Point2[0]]
                y += [Point1[1], Point2[1]]
                z += [Point1[2], Point2[2]]

        else:
            if len(RayListHistory[k]) > maxRays:
                rays_to_render = np.random.choice(RayListHistory[k], maxRays, replace=False)
            else:
                rays_to_render = RayListHistory[k]

            for j in rays_to_render:
                Point = j.point
                Vector = j.vector
                x += [Point[0], Point[0] + Vector[0] * EndDistance]
                y += [Point[1], Point[1] + Vector[1] * EndDistance]
                z += [Point[2], Point[2] + Vector[2] * EndDistance]
        points = np.column_stack((x, y, z))
        meshes += [pv.line_segments_from_points(points)]
    return meshes

def RayRenderGraph(OpticalChain, 
                   EndDistance=None, 
                   maxRays=300, 
                   OEpoints=2000, 
                   draw_mesh=False, 
                   cycle_ray_colors = False,
                   impact_points = False,
                   DrawDetectors=True,
                   DetectedRays = False,
                   Observers = dict()):
    """
    Renders an image of the Optical setup and the traced rays.

    Parameters
    ----------
        OpticalChain : OpticalChain
            List of objects of the ModuleOpticalOpticalChain.OpticalChain-class.

        EndDistance : float, optional
            The rays of the last ray bundle are drawn with a length given by EndDistance (in mm). If not specified,
            this distance is set to that between the source point and the 1st optical element.

        maxRays: int
            The maximum number of rays to render. Rendering all the traced rays is a insufferable resource hog
            and not required for a nice image. Default is 150.

        OEpoints : int
            How many little spheres to draw to represent the optical elements.  Default is 2000.

    Returns
    -------
        fig : Pyvista-figure-handle.
            Shows the figure.
    """

    RayListHistory = [OpticalChain.source_rays] + OpticalChain.get_output_rays()

    if EndDistance is None:
        EndDistance = np.linalg.norm(OpticalChain.source_rays[0].point - OpticalChain.optical_elements[0].position)

    print("...rendering image of optical chain...", end="", flush=True)
    fig = pvqt.BackgroundPlotter(window_size=(1500, 500), notebook=False) # Opening a window
    fig.set_background('white')
    
    if cycle_ray_colors:
        colors = mpu.generate_distinct_colors(len(OpticalChain)+1)
    else:
        colors = [[0.7, 0, 0]]*(len(OpticalChain)+1) # Default color: dark red

    # Optics display
    # For each optic we will send the figure to the function _RenderOpticalElement and it will add the optic to the figure
    for i,OE in enumerate(OpticalChain.optical_elements):
        color = pv.Color(colors[i+1])
        rgb = color.float_rgb
        h, l, s = rgb_to_hls(*rgb)
        s = max(0, min(1, s * 0.3))  # Decrease saturation
        l = max(0, min(1, l + 0.1))  # Increase lightness
        new_rgb = hls_to_rgb(h, l, s)
        darkened_color = pv.Color(new_rgb)
        mpm._RenderOpticalElement(fig, OE, OEpoints, draw_mesh, darkened_color, index=i)
    ray_meshes = _RenderRays(RayListHistory, EndDistance, maxRays)
    for i,ray in enumerate(ray_meshes):
        color = pv.Color(colors[i])
        fig.add_mesh(ray, color=color, name=f"RayBundle_{i}")
    if impact_points:
        for i,rays in enumerate(RayListHistory):
            points = np.array([list(r.point) for r in rays], dtype=np.float32)
            points = pv.PolyData(points)
            color = pv.Color(colors[i-1])
            fig.add_mesh(points, color=color, point_size=5, name=f"RayImpactPoints_{i}")
    
    detector_copies = {key: copy(OpticalChain.detectors[key][0]) for key in OpticalChain.detectors.keys()}
    detector_meshes_list = []
    detectedpoint_meshes = dict()
    
    if OpticalChain.detectors is not None and DrawDetectors:
        # Detector display
        for key in OpticalChain.detectors.keys():
            det = detector_copies[key]
            index = OpticalChain.detectors[key][1]
            if key in Observers:
                det.distance = Observers[key].value
                #Observers[key].register_calculation(lambda x: det.set_distance(x))
            mpm._RenderDetector(fig, det, name = key, detector_meshes = detector_meshes_list)
            if DetectedRays:
                RayListAnalysed = OpticalChain.get_output_rays()[index]
                points = det.get_3D_points(RayListAnalysed)
                points = pv.PolyData(points)
                detectedpoint_meshes[key] = points
                fig.add_mesh(points, color='purple', point_size=5, name=f"DetectedRays_{key}")
    detector_meshes = dict(zip(OpticalChain.detectors.keys(), detector_meshes_list))
    
    # Now we define a function that will move on the plot the detector with name "detname" when it's called
    def move_detector(detname, new_value):
        nonlocal fig, detector_meshes, detectedpoint_meshes, DetectedRays, detectedpoint_meshes, detector_copies, OpticalChain
        det = detector_copies[detname]
        index = OpticalChain.detectors[detname][1]
        det_mesh = detector_meshes[detname]
        translation = det.normal * (det.distance - new_value)
        det_mesh.translate(translation, inplace=True)
        det.distance = new_value
        if DetectedRays:
            points_mesh = detectedpoint_meshes[detname]
            points_mesh.points = det.get_3D_points(OpticalChain.get_output_rays()[index])
        fig.show()
    
    # Now we register the function to the observers
    for key in OpticalChain.detectors.keys():
        if key in Observers:
            Observers[key].register(lambda x: move_detector(key, x))

    #pv.save_meshio('optics.inp', pointcloud)  
    print(
        "\r\033[K", end="", flush=True
    )  # move to beginning of the line with \r and then delete the whole line with \033[K
    fig.show()
    return fig

moc.OpticalChain.render = RayRenderGraph

# %% Asphericity

def plot_asphericity(Mirror, Npoints=1000):
    """
    This function displays a map of the asphericity of the mirror.
    It's a scatter plot of the points of the mirror surface, with the color representing the distance to the closest sphere.
    The closest sphere is calculated by the function get_closest_sphere, so least square method.

    Parameters
    ----------
    Mirror : Mirror
        The mirror to analyse.

    Npoints : int, optional
        The number of points to sample on the mirror surface. The default is 1000.
    
    Returns
    -------
    fig : Figure
        The figure of the plot.
    """
    plt.ion()
    fig = plt.figure()
    ax = Mirror.support._ContourSupport(fig)
    center, radius = man.get_closest_sphere(Mirror, Npoints)
    Points = mpm.sample_support(Mirror.support, Npoints=1000)
    Points += Mirror.r0[:2]
    Z = Mirror._zfunc(Points)
    Points = mgeo.PointArray([Points[:, 0], Points[:, 1], Z]).T
    X, Y = Points[:, 0] - Mirror.r0[0], Points[:, 1] - Mirror.r0[1]
    Points_centered = Points - center
    Distance = np.linalg.norm(Points_centered, axis=1) - radius
    Distance*=1e3 # To convert to µm
    p = plt.scatter(X, Y, c=Distance, s=15)
    divider = man.make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(p, cax=cax)
    cbar.set_label("Distance to closest sphere (µm)")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.title("Asphericity map", loc="right")
    plt.tight_layout()

    bbox = ax.get_position()
    bbox.set_points(bbox.get_points() - np.array([[0.01, 0], [0.01, 0]]))
    ax.set_position(bbox)
    plt.show()
    return fig


# %% Caustics

def plot_caustics(OpticalChain, Range, DetectorName="Focus" , Npoints=1000, Nrays=1000):
    """
    This function displays the caustics of the rays on the detector.
    To do so, it calculates the intersections of the rays with the detector over a 
    range determined by the parameter Range, and then plots the standard deviation of the
    positions in the x and y directions.

    Parameters
    ----------
    OpticalChain : OpticalChain
        The optical chain to analyse.

    DetectorName : str
        The name of the detector on which the caustics are calculated.
    
    Range : float
        The range of the detector over which to calculate the caustics.

    Npoints : int, optional
        The number of points to sample on the detector. The default is 1000.
    
    Returns
    -------
    fig : Figure
        The figure of the plot.
    """
    distances = np.linspace(-Range, Range, Npoints)
    Det,Index = OpticalChain.detectors[DetectorName]
    Rays = OpticalChain.get_output_rays()[Index]
    Rays = np.random.choice(Rays, Nrays, replace=False)
    LocalRayList = [r.to_basis(*Det.basis) for r in Rays]
    Points = mgeo.IntersectionRayListZPlane(LocalRayList, distances)
    x_std = []
    y_std = []
    for i in range(len(distances)):
        x_std.append(mp.StandardDeviation(Points[i][:,0]))
        y_std.append(mp.StandardDeviation(Points[i][:,1]))
    plt.ion()
    fig, ax = plt.subplots()
    ax.plot(distances, x_std, label="x std")
    ax.plot(distances, y_std, label="y std")
    ax.set_xlabel("Detector distance (mm)")
    ax.set_ylabel("Standard deviation (mm)")
    ax.legend()
    plt.title("Caustics")
    plt.show()
    return fig