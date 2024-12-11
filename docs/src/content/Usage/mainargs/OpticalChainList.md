---
title: "OpticalChainList"
weight: 1
type: docs
---

{{% alert title="Summary" color="primary" %}}
Use the `OEPlacement()` function on a list of optical elements, a `SourceProperties` dictionnary, a list of distances and a list of incidence angles.
{{% /alert %}}

## Explanation
The first argument of `main()` is `OpticalChainList`, a complete description of the optical setup that needs to be simulated. This includes the absolute positions and orientations of all the optical elements apart from the detector.

Positionning these elements in 3D space by hand would be unnecessarily complicated. Instead, a utility function `OEPlacement()` has been written to allow you to automatically position these elements, provided you specify the distances between them and the incidence angle on each element.

To position these elements, a prealignment is performed, much like with an alignment laser in real life. A single ray (the central ray within the actual light source) is cast and the Optical Elements are successively positionned in the path of that ray at the specified distance and incidence angle. A reflected ray is then calculated and the algorithm uses it to add the next Optical Element.

{{% alert title="FYI" color="info" %}}
Using mirrors with deformations would previously skew the algorithm as the central ray would get randomly reflected by the deformation in the center of the mirror. This is no longer and issue as during prealignment in the `OEPlacement()` function, the mirror defects are automatically turned off.
{{% /alert %}}

It's perfectly possible to modify the alignment after this prealignment procedure. See for instance [here](/examples/config_2toroidals_f-x-f/).

## Using `OEPlacement()`

To use `OEPlacement()` we logically need:
 - A Python list of mirrors or masks to be positionned
 - A source of light to be used for prealignment (a virtual HeNe laser) described by a `dict` (see [SourceProperties](/usage/mainargs/sourceproperties/))
 - A list of distances between optics. The first value is the distance between source and first optic.
 - A list of incidence angles of the central ray on the optical elements.
