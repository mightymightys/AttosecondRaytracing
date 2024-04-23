"""
Provides functions for analysis of the output ray-bundles calculated through ray-tracing, and also for the ART's standard visualization options.

These are called by ARTmain according to the values in the *AnalysisOptions*-dictionary filled in the CONFIG-scripts.




Created in Apr 2020

@author: Anthony Guillaume + Stefan Haessler
"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
import pyvistaqt as pvqt
import colorcet as cc
import colorsys

import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleGeometry as mgeo
import itertools


# %%
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

    DetectorPointList2DCentre = Detector.get_PointList2DCentre(RayListAnalysed)
    DectectorPoint2D_Xcoord = np.array([k[0] * 1e3 for k in DetectorPointList2DCentre])  # mm to µm
    DectectorPoint2D_Ycoord = np.array([k[1] * 1e3 for k in DetectorPointList2DCentre])

    FocalSpotSize = mgeo.DiameterPointList(DetectorPointList2DCentre)
    SpotSizeSD = mp.StandardDeviation(DetectorPointList2DCentre)
    return DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD


# %%
def getETransmission(RayListIn, RayListOut) -> float:
    """
    Calculates the energy transmission from RayListIn to RayListOut in percent by summing up the
    intensity-property of the individual rays.

    Parameters
    ----------
        RayListIn, RayListOut : list(Ray)
            Lists of objects of the ModuleOpticalRay.Ray-class.

    Returns
    -------
        ETransmission : float
    """
    ETransmission = 100 * sum(Ray.intensity for Ray in RayListOut) / sum(Ray.intensity for Ray in RayListIn)
    return ETransmission


# %%
def GetResultSummary(Detector, RayListAnalysed, verbose=False):
    """
    Calculate and return FocalSpotSize-standard-deviation and Duration-standard-deviation
    for the given Detector and RayList.
    If verbose, then also print a summary of the results for the given Detector.

    Parameters
    ----------
        Detector : Detector
            An object of the ModuleDetector.Detector-class.

        RayListAnalysed : list(Ray)
            List of objects of the ModuleOpticalRay.Ray-class.

        verbose : bool
            Whether to print a result summary.

    Returns
    -------
        FocalSpotSizeSD : float

        DurationSD : float
    """
    DetectorPointList2DCentre = Detector.get_PointList2DCentre(RayListAnalysed)
    FocalSpotSizeSD = mp.StandardDeviation(DetectorPointList2DCentre)
    DelayList = Detector.get_Delays(RayListAnalysed)
    DurationSD = mp.StandardDeviation(DelayList)

    if verbose:
        FocalSpotSize = mgeo.DiameterPointList(DetectorPointList2DCentre)
        summarystring = (
            "At the detector distance of "
            + "{:.3f}".format(Detector.get_distance())
            + " mm we get:\n"
            + "Spatial std : "
            + "{:.3f}".format(FocalSpotSizeSD * 1e3)
            + " \u03BCm and min-max: "
            + "{:.3f}".format(FocalSpotSize * 1e3)
            + " \u03BCm\n"
            + "Temporal std : "
            + "{:.3e}".format(DurationSD)
            + " fs and min-max : "
            + "{:.3e}".format(max(DelayList) - min(DelayList))
            + " fs"
        )

        print(summarystring)

    return FocalSpotSizeSD, DurationSD


# %%
def SpotDiagram(RayListAnalysed, Detector, DrawAiryAndFourier=False, ColorCoded=None) -> plt.Figure:
    """
    Produce a an interactive figure with the spot diagram resulting from the RayListAnalysed
    hitting the Detector.
    The detector distance can be shifted with the left-right cursor keys.
    If DrawAiryAndFourier is True, a circle with the Airy-spot-size will be shown.
    The 'spots' can optionally be color-coded by specifying ColorCoded
    as ["Intensity","Incidence","Delay"].

    Parameters
    ----------
        RayListAnalysed : list(Ray)
            List of objects of the ModuleOpticalRay.Ray-class.

        Detector : Detector
            An object of the ModuleDetector.Detector-class.

        DrawAiryAndFourier : bool, optional
            Whether to draw a circle with the Airy-spot-size. The default is False.

        ColorCoded : str, optional
            Color-code the spots according to one of ["Intensity","Incidence","Delay"]. The default is None.

    Returns
    -------
        fig : matlplotlib-figure-handle.
            Shows the interactive figure.
    """
    NumericalAperture = mp.ReturnNumericalAperture(RayListAnalysed, 1)  # NA determined from final ray bundle
    Wavelength = RayListAnalysed[0].wavelength
    if DrawAiryAndFourier:
        AiryRadius = mp.ReturnAiryRadius(Wavelength, NumericalAperture) * 1e3  # in µm
    else:
        AiryRadius = 0

    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = _getDetectorPoints(
        RayListAnalysed, Detector
    )

    if ColorCoded == "Intensity":
        IntensityList = [k.intensity for k in RayListAnalysed]
        z = np.asarray(IntensityList)
        zlabel = "Intensity (arb.u.)"
        title = "Intensity + Spot Diagram\n press left/right to move detector position"
        addLine = ""
    elif ColorCoded == "Incidence":
        IncidenceList = [np.rad2deg(k.incidence) for k in RayListAnalysed]  # degree
        z = np.asarray(IncidenceList)
        zlabel = "Incidence angle (deg)"
        title = "Ray Incidence + Spot Diagram\n press left/right to move detector position"
        addLine = ""
    elif ColorCoded == "Delay":
        DelayList = Detector.get_Delays(RayListAnalysed)
        DurationSD = mp.StandardDeviation(DelayList)
        z = np.asarray(DelayList)
        zlabel = "Delay (fs)"
        title = "Delay + Spot Diagram\n press left/right to move detector position"
        addLine = "\n" + "{:.2f}".format(DurationSD) + " fs SD"
    else:
        z = "red"
        title = "Spot Diagram\n press left/right to move detector position"
        addLine = ""

    Dist = Detector.get_distance()
    distStep = min(50, max(0.0005, round(FocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000))  # in mm

    plt.ion()
    fig, ax = plt.subplots()
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
        label="{:.3f}".format(Dist) + " mm\n" + "{:.1f}".format(SpotSizeSD * 1e3) + " \u03BCm SD" + addLine,
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

    movingDetector = Detector.copy_detector()

    def press(event):
        nonlocal Dist, distStep, movingDetector, ColorCoded, zlabel, cbar
        if event.key == "right":
            movingDetector.shiftByDistance(distStep)
            Dist += distStep
        elif event.key == "left":
            if Dist > 1.5 * distStep:
                movingDetector.shiftByDistance(-distStep)
                Dist -= distStep
            else:
                movingDetector.shiftToDistance(0.5 * distStep)
                Dist = 0.5 * distStep
        else:
            return None
        newDectectorPoint2D_Xcoord, newDectectorPoint2D_Ycoord, newFocalSpotSize, newSpotSizeSD = _getDetectorPoints(
            RayListAnalysed, movingDetector
        )
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

        foo.set_label(
            "{:.3f}".format(Dist) + " mm\n" + "{:.1f}".format(newSpotSizeSD * 1e3) + " \u03BCm SD" + newaddLine
        )
        ax.legend(loc="upper right")

        axisLim = 1.1 * max(AiryRadius, 0.5 * newFocalSpotSize * 1000)
        ax.set_xlim(-axisLim, axisLim)
        ax.set_ylim(-axisLim, axisLim)

        distStep = min(
            50, max(0.0005, round(newFocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000)
        )  # in mm

        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("key_press_event", press)

    plt.show()

    return fig


# %%
def _drawDelayGraph(RayListAnalysed, Detector, Distance, DeltaFT, DrawAiryAndFourier=False, ColorCoded=None, fig=None):
    """
    Draws the 3D-delay-spot-diagram for a fixed detector. See more doc in the function below.
    """
    NumericalAperture = mp.ReturnNumericalAperture(RayListAnalysed, 1)  # NA determined from final ray bundle
    Wavelength = RayListAnalysed[0].wavelength
    AiryRadius = mp.ReturnAiryRadius(Wavelength, NumericalAperture) * 1e3  # in µm

    DectectorPoint2D_Xcoord, DectectorPoint2D_Ycoord, FocalSpotSize, SpotSizeSD = _getDetectorPoints(
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


# %%
def DelayGraph(
    RayListAnalysed, Detector, DeltaFT: (int, float), DrawAiryAndFourier=False, ColorCoded=None
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
    Dist = Detector.get_distance()
    fig, NumericalAperture, AiryRadius, FocalSpotSize = _drawDelayGraph(
        RayListAnalysed, Detector, Dist, DeltaFT, DrawAiryAndFourier, ColorCoded
    )

    distStep = min(50, max(0.0005, round(FocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000))  # in mm

    movingDetector = Detector.copy_detector()

    def press(event):
        nonlocal Dist, distStep, movingDetector, fig
        if event.key == "right":
            movingDetector.shiftByDistance(distStep)
            Dist += distStep
            ax = fig.axes[0]
            cam = [ax.azim, ax.elev, ax._dist]
            fig, sameNumericalAperture, sameAiryRadius, newFocalSpotSize = _drawDelayGraph(
                RayListAnalysed, movingDetector, Dist, DeltaFT, DrawAiryAndFourier, ColorCoded, fig
            )
            ax = fig.axes[0]
            ax.azim, ax.elev, ax._dist = cam
        elif event.key == "left":
            if Dist > 1.5 * distStep:
                movingDetector.shiftByDistance(-distStep)
                Dist -= distStep
            else:
                movingDetector.shiftToDistance(0.5 * distStep)
                Dist = 0.5 * distStep
            ax = fig.axes[0]
            cam = [ax.azim, ax.elev, ax._dist]

            fig, sameNumericalAperture, sameAiryRadius, newFocalSpotSize = _drawDelayGraph(
                RayListAnalysed, movingDetector, Dist, DeltaFT, DrawAiryAndFourier, ColorCoded, fig
            )
            ax = fig.axes[0]
            ax.azim, ax.elev, ax._dist = cam
        else:
            return fig
        distStep = min(
            50, max(0.0005, round(newFocalSpotSize / 8 / np.arcsin(NumericalAperture) * 10000) / 10000)
        )  # in mm

    fig.canvas.mpl_connect("key_press_event", press)

    return fig


# %%
def MirrorProjection(OpticalChain, ReflectionNumber: int, Detector=None, ColorCoded=None) -> plt.Figure:
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

    Position = OpticalChain.optical_elements[ReflectionNumber].position
    n = OpticalChain.optical_elements[ReflectionNumber].normal
    m = OpticalChain.optical_elements[ReflectionNumber].majoraxis

    RayListAnalysed = OpticalChain.get_output_rays()[ReflectionNumber]
    # transform rays into the mirror-support reference frame
    # (same as mirror frame but without the shift by mirror-centre)
    RayList = mgeo.TranslationRayList(RayListAnalysed, -Position)
    RayList = mgeo.RotationRayList(RayList, n, np.array([0, 0, 1]))
    mPrime = mgeo.RotationPoint(m, n, np.array([0, 0, 1]))
    RayList = mgeo.RotationRayList(RayList, mPrime, np.array([1, 0, 0]))

    x = np.asarray([k.point[0] for k in RayList])
    y = np.asarray([k.point[1] for k in RayList])
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
        if Detector is not None:
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
    ax = OpticalChain.optical_elements[ReflectionNumber].type.support._ContourSupport(fig)
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


# %%
def _RenderOpticalElement(OE, OEpoints, draw_mesh = False):
    OpticPointList, edge_faces = OE.type.get_grid3D(OEpoints, edges=True)  # in the optic's coordinate system
    # transform OpticPointList into "lab-frame"
    OpticPointList = mgeo.TranslationPointList(OpticPointList, -OE.type.get_centre())
    MirrorMajorAxisPrime = mgeo.RotationPoint(OE.majoraxis, OE.normal, np.array([0, 0, 1]))
    OpticPointList = mgeo.RotationPointList(OpticPointList, np.array([1, 0, 0]), MirrorMajorAxisPrime)
    OpticPointList = mgeo.RotationPointList(OpticPointList, np.array([0, 0, 1]), OE.normal)
    OpticPointList = mgeo.TranslationPointList(OpticPointList, OE.position)
    
    # shift optics surfaces a little back normally?
    if not draw_mesh:    
        OpticPointList = [point - OE.normal*0.5 for point in OpticPointList]

    optic_pts = pv.PolyData(OpticPointList)
    tess = None
    if draw_mesh:
        pts_coord = pv.PolyData(OpticPointList)
        lines = list(
            itertools.chain.from_iterable([[[2, e[i], e[i + 1]] for i in range(len(e) - 1)] for e in edge_faces])
        )
        faces = list(itertools.chain.from_iterable([[len(i) - 1] + i[:-1] for i in edge_faces]))
        if lines == []:
            lines = [0]
        if faces == []:
            faces = [0]
        edges = pv.PolyData(
            OpticPointList,
            lines=lines,
            faces=faces,
        )
        # Can't figure out some edge cases such as when part of the support is outside of the mirror
        tess = pts_coord.delaunay_2d(edge_source=edges)
    return optic_pts, tess

def _RenderRays(RayListHistory, EndDistance, maxRays=150, color_by_number = True):
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

def generate_distinct_colors(num_colors):
    # Get a color palette from colorcet
    palette = cc.glasbey

    # Make sure the number of colors does not exceed the palette length
    num_colors = min(num_colors, len(palette))

    # Slice the palette to get the desired number of colors
    distinct_colors = palette[:num_colors]

    return distinct_colors

def RayRenderGraph(OpticalChain, EndDistance=None, maxRays=300, OEpoints=3000, scale_spheres=5.0, draw_mesh=False, cycle_ray_colors = False):
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
    fig = pvqt.BackgroundPlotter(window_size=(1500, 500), notebook=False)
    fig.set_background('white')
    
    ray_meshes = _RenderRays(RayListHistory, EndDistance, maxRays)
    if cycle_ray_colors:
        colors = generate_distinct_colors(len(ray_meshes))
    else:
        colors = [[0.7, 0, 0]]*len(ray_meshes)
    for i,ray in enumerate(ray_meshes):
        color = pv.Color(colors[i])
        fig.add_mesh(ray, color=color)

    # Optics display
    for i,OE in enumerate(OpticalChain.optical_elements):
        pointcloud, mesh = _RenderOpticalElement(OE, OEpoints, draw_mesh)
        color = pv.Color(colors[i+1])
        color = colorsys.hsv_to_rgb(*(colorsys.rgb_to_hsv(*color[:-1]))*np.array([1,0.2,1]))
        if draw_mesh and mesh is not None:
            fig.add_mesh(mesh, color = color)
        fig.add_mesh(pointcloud, point_size=scale_spheres, color = color, render_points_as_spheres = True)
    fig.show()
    
    print(
        "\r\033[K", end="", flush=True
    )  # move to beginning of the line with \r and then delete the whole line with \033[K
    return fig




def show():
    plt.show(block=False)
