---
title: "AnalysisOptions"
weight: 4
type: docs
---

{{% alert title="Summary" color="primary" %}}
`AnalysisOptions` is a python dictionnary
{{% /alert %}}

```Python
DefaultAnalysisOptions = {
    "verbose": True,  # print intermediate results and info in the console?
    "plot_Render": False,  # render optical elements and rays, and how many rays to render?
    "maxRaysToRender": 150,
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
    "OEPointsToRender": 2000,
    "OEPointsScale": 0.5,
}
```
