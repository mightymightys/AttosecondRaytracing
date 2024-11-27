This package contains the core code required for the ART - AttosecondRaytracing code. It contains the definitions of basic optical elements as well as the actual ray-tracing code. If you want to develop your own visualisation/GUI for ART, this is the package to install.
In normal use case, it's installed as a dependency of [AttosecondRayTracing](https://pypi.org/project/AttosecondRayTracing/) which is the user facing code.

[ART - Attosecond Ray Tracing](https://github.com/mightymightys/AttosecondRaytracing) - is a free python code written by Stefan Haessler, André Kalouguine, and Anthony Guillaume of
[Laboratoire d'Optique Appliquée (LOA), CNRS, Institut Polytechnique de Paris, France](https://loa.ensta-paris.fr/research/pco-research-group/)
and Charles Bourassin-Bouchet [Laboratoire Charles Fabry (LCF)), Institut d’Optique, CNRS, Université Paris-Saclay, France](https://www.lcf.institutoptique.fr/en/groups/optique-xuv).

It does ray tracing calculations with a specific focus on ultrashort (femto- and attosecond) laser pulses.
Therefore the code currently focuses on reflective optics, freely arrangeable including grazing incidence configurations.

Ten years ago, Charles made geometric optics calculations that demonstrated how sensitive attosecond pulses
are to spatio-temporal distortions and how easily such distortions are picked up in the reflective
grazing-incidence optical setups required to transport and refocus them
[[C. Bourassin-Bouchet et al. “How to focus an attosecond pulse”. Opt.
Express 21, 2506 (2013)](http://dx.doi.org/10.1364/oe.21.002506); [C. Bourassin-Bouchet et al. “Spatiotemporal distortions of
attosecond pulses”. JOSA A 27, 1395 (2010)](https://www.osapublishing.org/josaa/abstract.cfm?uri=josaa-27-6-1395)].
ART now makes these calculations avaible to the ultrafast optics community in a free and (hopefully) easily accessible python code.

A publication of simulations of the beam transport and focusing of high-numerical-aperture XUV attosecond pulses
is in preparation and will become the reference that we ask you to cite if you have used this code in your work.

The detailed documentaion can be found [here](https://mightymightys.github.io/AttosecondRaytracing/).
