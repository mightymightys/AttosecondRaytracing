---
title: "Installation"
linkTitle: "Installation"
type: docs
weight: 2
---

The code is split into two Python packages:
 - `AttosecondRayTracing-core` that provides the `ARTcore` module. It contains the actual ray-tracing code, the definitions of the mirrors, sources and detectors.
 - `AttosecondRayTracing` that provides the `ART` module and requires `AttosecondRayTracing-core` as a dependency. It contains the plotting and visualisation functions as well as the analsysis functions. It's the most user-facing part of the codebase. 

Unless you really only need the functionality of `AttosecondRayTracing-core`, we recommend installing `AttosecondRayTracing` that will pull in the required dependencies.

## Basic installation
### Using `pip`
If installing using `pip`, we recommend installing the dependencies in a virtual environment, for instance using 
```Shell
python -m venv <new_virtual_environment_folder>
```
This lets you install and use the software without interfering with the system installation of Python.

The installation of the package itself couldn't be easier:
```Shell
pip install AttosecondRayTracing
```

### Using Anaconda
Just as with `pip`, we recommend using a separate virtual environment to install and use ART. 
TODO actual installation instructions.

## Developper installation
If you want to contribute to the development of ART or simply want to modify something for your own needs, you can directly download the source code from the github page.

For instance:
```Shell
git clone https://github.com/mightymightys/AttosecondRaytracing.git
```

### Contributing
We appreciate any contributions to the code. The best way to do so is by making pull requests on Github.