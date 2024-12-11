---
title: "CONFIG_2toroidals_f-x-f.py"
weight: 2
type: docs
---

This example configuration describes a refocusing setup comprised of two identical toroidal mirrors with a focal length of 500 and an angle of incidence of 80°. The beam from a point source is first cut to 35mrad by a circular mask, it then impacts the first toroidal mirror positionned one focal length away from the source. A second toroidal mirror is positionned at a variable distance away from the first one, with

![Refocusing using 2 toroidal mirrors](../CONFIG_2toroidals_f-x-f.svg)

You can download the configuration file here: [``CONFIG_2toroidals_f-x-f.py``](../examples/CONFIG_2toroidals_f-x-f.py)

Opening the configuration code, the first lines import the necessary modules, mostly from `ARTcore`.
```Python
import numpy as np
import ARTcore.ModuleMirror as mmirror
import ARTcore.ModuleSupport as msupp
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleMask as mmask
from ART.ARTmain import main
```

The actual raytracing calculations are performed by the [`RayTracingCalculation`](artcoreapi/ModuleProcessing.html#RayTracingCalculation) function. It requires a list of source rays and a list of optics placed in 3D space and returns a list of output rays.

However, positionning the elements in 3D space is cumbersome, so instead we will use the [`OEPlacement()`](/artcoreapi/ModuleProcessing.html#OEPlacement) function. It will automatically position the optical elements in 3D space in a way that respects the distances and incidence angles that the user specifies. It outputs a container object for the entire setup: an  [`OpticalChain`](artcoreapi/ModuleOpticalChain.html#OpticalChain).

To position the elements, the `OEPlacement()` function itself performs a raytracing calculation. Instead of accepting directly a source of rays, the function takes in a `dict` describing the properties of the source and generates its own list of rays

            

### Source
The source is nothing more than a list of [`Ray`](artcoreapi/ModuleOpticalRay.html#Ray). The [`ModuleSource`](artcoreapi/ModuleSource) module provides some utility functions that generate some useful examples of sources. In this case, we will use a simple point source with the following properties:
```Python
SourceProperties = {
    'Divergence' : 50e-3/2, # half-angle in rad
    'SourceSize' : 0, # diameter in mm
    'Wavelength' : 80e-6, # 40 nm
    'DeltaFT'    : 0.5,   # in fs
    'NumberRays' : 1000   # launch 1000 rays in the beginning
}
```
