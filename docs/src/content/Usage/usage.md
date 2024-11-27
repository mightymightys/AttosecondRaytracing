---
title: "Usage"
linkTitle: "Usage"
type: docs
weight: 2
---

Before starting, we strongly recommend checking out the [conventions](../conventions) we use in ART. 

The entry point into the program is the `main()` function found in the `ARTmain.py` file. It has five arguments, one of which is optional:
 - `OpticalChainList`: An `OpticalChainList` object containing all the optical elements to be simulated with their positions as well as a list of initial source rays.
 - `SourceProperties`: A `dict` of the properties of the source that will be used for the light simulation
 - `DetectorOptions`: A `dict` describing the position of the light ray detector or the parameters allowing it to be auto-positionned.
 - `AnalysisOptions`: A `dict` specifying which plots to display and their options.
 - `save_file_name` (optional): A file name for saving the raw data to be analysed later.

Thus, running the program requires constructing these objects and calling `main()`. There are two ways to do so:
 - Construct them in a script or in an interactive shell or IPython notebook and call the function.
 - Run the `ARTmain.py` file as a program and pass a config file as its first argument.


## Using a configuration file

In practice, the config file is simply a Python script similar to what one would write using the first way. The main difference is that it doesn't contain a call to the `main()` function. 

Indeed, running `ARTmain.py` as a standalone program with a config file simply loads the file as a module and checks for the presence of attributes with the same names as the arguments of `main()`.

Thus, to write a config file, you simply have to create variables with the same names as the arguments of `main()`.

## Running directly

Same as with using a config file, you need to create the arguments for the `main()` function. Calling it will run the program as well as return the raw data for further analysis.

