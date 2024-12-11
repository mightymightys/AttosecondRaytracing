![A rendering of two toroidal mirrors with an intermediate collimated section.](./docs/pdoc/images/doc_illustrationrender.png)


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

A quickstart documentation will very soon (promised ;-) be added in the following.

## Dependencies

The code requires Python 3.6 or newer and depends on the libraries [NumPy](https://numpy.org), 
[Numpy-Quaternion](https://github.com/moble/quaternion),  [matplotlib](https://matplotlib.org),
and for 3D-rendering of optical configurations and rays,  [PyVista](https://github.com/pyvista/pyvista).

### Using `pip`
If installing using `pip`, we recommend installing the dependencies in a virtual environment, for instance using 
```Shell
python -m venv <new_virtual_environment_folder>
```
This lets you install and use the software without interfering with the system installation of Python.

As for the installation of the dependencies:
```Shell
pip install numpy numpy-quaternion matplotlib pyvista scipy colorcet
```

### Using Anaconda
Just as with `pip`, we recommend using a separate virtual environment to install and use ART. The dependencies are all available in the *conda-forge* repository.

## Installation

You can just download the code as a zip file from here. Or if you are using the git version control software,
ou can clone the repository like so:

```Shell
git clone https://github.com/mightymightys/AttosecondRaytracing.git
```

You are welcome to fork the code and contribute to its further development


## Running ART 

To run ART, you run the appropriately named ***ARTmain.py*** in the console, supplying a
configuration-file as an argument, like so:
  
    python ARTmain.py CONFIG_xxxx.py

The configuration file is itself written in python, and the user is responsible for enuring
that it doesn’t contain any harmful code. A template for such a configuration file is given
by ***CONFIG_template.py*** and a number of example configurations files are provided to
be tried, explored and adapted.

Alternatively, the user may use the configuration file as a launch script, and run the
main program as indicated at the end of the template and example configuration files.
This is practicle, e.g., when writing the configuration file in an IDE like *Spyder*, which
features an IPython-console and let the user run the configuration directly.