"""
![A rendering of two toroidal mirrors with an intermediate collimated section.](../docs/doc_illustrationrender.png)
                                                          
ART - Attosecond Ray Tracing - is a free python code written by Stefan Haessler, Anthony Guillaume of
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

### Installation / Dependencies

You can just download the code as a zip file from here. Or if you are using the git version control software,
ou can clone the repository like so:
    
    git clone https://github.com/mightymightys/AttosecondRaytracing.git

You are welcome to fork the code and contribute to its further development!

The code requires Python 3.6 or newer and depends on the libraries [NumPy](https://numpy.org), 
[Numpy-Quaternion](https://github.com/moble/quaternion),  [matplotlib](https://matplotlib.org),
and for 3D-rendering of optical configurations and rays,  [Mayavi](https://docs.enthought.com/mayavi/mayavi).

We strongly recommend using an anaconda/miniconda python distribution, and recreating the virtual
conda-environment ***ARTenv***, fixed in the file *ARTenvironment.yml* contained in the repository.
This will make sure you have a combination of versions of all dependencies that has been tested to work
as expected. In a terminal, do:
    
    conda env create -f ARTenvironment.yml
    conda activate ARTenv

Otherwise you can install the crucial dependencies easily if you use a miniconda distribution by entering in the Anaconda prompt:
    
    conda install -c anaconda numpy
    conda install -c anaconda matplotlib
    conda install -c conda-forge quaternion
    conda install -c conda-forge mayavi

With some luck, this will let ART work in the base environment. Howver, in particular the Mayavi package
(or rather the VTK package that it depends on) is not always available for the most recent python version.


### Running ART 

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

### More detailed instructions to come...

The [Optical Chain](ART/ModuleOpticalChain.html) is a class that represents the whole optical setup to be simulated,
and is the central object the user will deal with...

"""


