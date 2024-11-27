---
title: Attosecond Ray Tracing
linkTitle: Attosecond Ray Tracing
type: "docs"
---

![A rendering of two toroidal mirrors with an intermediate collimated section.](doc_illustrationrender.png)

{{% pageinfo %}}
This documentation web page is a work in progress. Please do not hesitate to raise issues in the [Github repository](https://github.com/mightymightys/AttosecondRaytracing) of the project.
{{% /pageinfo %}}

The AttosecondRayTracing code is a tool for the attosecond-physics community to assist with the design of beam transport and refocusing setups. 

It is written by Stefan Haessler, André Kalouguine, and Anthony Guillaume of
[Laboratoire d'Optique Appliquée (LOA), CNRS, Institut Polytechnique de Paris, France](https://loa.ensta-paris.fr/research/pco-research-group/)
and Charles Bourassin-Bouchet [Laboratoire Charles Fabry (LCF)), Institut d’Optique, CNRS, Université Paris-Saclay, France](https://www.lcf.institutoptique.fr/en/groups/optique-xuv).

It does ray tracing calculations for (a succession of) reflective optics, freely arrangeable including grazing incidence configurations.
It keeps track of the optical paths of each ray and therefore supplies precise information about the delays with which
the rays hit a virtual detector plane. These can --within the limits of geometric optics-- describe spatio-temporal distortions
introduced by aberrations of the optical setup.

Such calculations were already published by Charles ten years ago for single-mirror setups, demonstrating how
sensitive attosecond pulses are to spatio-temporal distortions and how easily such distortions are picked up in
the reflective grazing-incidence optical setups required to transport and refocus them
[[C. Bourassin-Bouchet et al. “How to focus an attosecond pulse”. Opt.
Express 21, 2506 (2013)](http://dx.doi.org/10.1364/oe.21.002506); [C. Bourassin-Bouchet et al. “Spatiotemporal distortions of
attosecond pulses”. JOSA A 27, 1395 (2010)](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-27-6-1395)].

A publication of simulations of the beam transport and focusing of high-numerical-aperture XUV attosecond pulses
with multi-mirror setups is in preparation and will become the reference that we ask you to cite if you have used
this code in your work.


Please start by reading the [Installation and Usage](usage) pages.
