---
title: "SourceProperties"
weight: 2
type: docs
---

{{% alert title="Summary" color="primary" %}}
`SourceProperties` is a python dictionnary
{{% /alert %}}

Running the simulation requires describing the light source. The properties that need to be defined are the following:
 - `"DeltaFT"`: This is not actually a property of the light source. Indeed, in the simulation, all photons are considered to be emitted simultaneously in an infinitely short burst. This parameter is actually an analysis parameter, used to check whether the light rays after the optical system have a short enough temporal spread.
 - `"Wavelength"`: Similarly to `DeltaFT`, this is actually an analysis parameter that is used as a reference for the spatial spread of the light rays after the optical system. It can however be used in the future to make wavelength-dependent optics such as gratings.
 - `"Divergence"`: This is an actual propertany of the light source: it's the half-angle of the cone of light emitted by the source as illustrated below.
 - `"SourceSize"`: Again, a property of the light source. It's the diameter of the light source (considered to be a disc). If the light source is a point source, it should be set to 0.
 - `"NumberRays"`: The number of rays that the light source will cast and that will be used for the simulation.

![Illustration of divergence of light source](/api/pointsource.svg)


Below we provide an example of a 5mm wide source of 80nm light with a 50mrad total divergence, casting 1000 rays and compared in the end to an Airy cylinder of 0.5 femtoseconds.
```Python
SourceProperties = {
    'Divergence' : 50e-3/2,
    'SourceSize' : 5,
    'Wavelength' : 80e-6,
    'DeltaFT'    : 0.5,
    'NumberRays' : 1000
}
```
