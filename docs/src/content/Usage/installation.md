---
title: "Installation"
linkTitle: "Installation"
weight: 2
---

## Dependencies
The code requires Python 3.6 or newer and depends on the libraries [NumPy](https://numpy.org), 
[Numpy-Quaternion](https://github.com/moble/quaternion),  [matplotlib](https://matplotlib.org),
and for 3D-rendering of optical configurations and rays,  [PyVista](https://github.com/pyvista/pyvista).
It also now requires SciPy ??

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
Just as with `pip`, we recommend using a separate virtual environment to install and use ART. 

## Installation

You can just download the code as a zip file from here. Or if you are using the git version control software,
ou can clone the repository like so:

```Shell
git clone https://github.com/mightymightys/AttosecondRaytracing.git
```

You are welcome to fork the code and contribute to its further development


