# -*- coding: utf-8 -*-
"""
Created in Jan 2023

@author: Stefan Haessler
"""
import os
import sys
from datetime import datetime
import numpy as np
import importlib.util
import ARTcore.ModuleProcessing as mp
import ARTcore.ModuleDetector as mdet
import ART.ModuleAnalysisAndPlots as mplots
#import matplotlib.pyplot as plt
import ARTcore.ModuleOpticalChain as moc


def print_banner(i):
    """Show the important name banner and the version-number and date."""

    banner = [
        r"""
    _  _   _                              _    ___             _____            _
   /_\| |_| |_ ___ ___ ___ __ ___ _ _  __| |  | _ \__ _ _  _  |_   _| _ __ _ __(_)_ _  __ _
  / _ \  _|  _/ _ (_-</ -_) _/ _ \ ' \/ _` |  |   / _` | || |   | || '_/ _` / _| | ' \/ _` |
 /_/ \_\__|\__\___/__/\___\__\___/_||_\__,_|  |_|_\__,_|\_, |   |_||_| \__,_\__|_|_||_\__, |
                                                        |__/                          |___/"""
    ]

    banner.append(
        r"""
   ___  __  __                                __   ___              ______             _
  / _ |/ /_/ /____  ___ ___ _______  ___  ___/ /  / _ \___  __ __  /_  __/______  ____(_)__  ___ _
 / __ / __/ __/ _ \(_-</ -_) __/ _ \/ _ \/ _  /  / , _/ _ `/ // /   / / / __/ _ `/ __/ / _ \/ _ `/
/_/ |_\__/\__/\___/___/\__/\__/\___/_//_/\_,_/  /_/|_|\_,_/\_, /   /_/ /_/  \_,_/\__/_/_//_/\_, /
                                                          /___/                            /___/  """
    )

    niceline = "___________________________________________________________________________________________________"

    if abs(i) > len(banner) - 1:
        i = -1
    print(niceline)
    print(banner[i], flush=True)

    with open("ART/VERSION", "r") as fvers:
        __version__ = fvers.read().strip()

    mod_time = os.path.getmtime("ART/VERSION")
    mod_time_str = datetime.fromtimestamp(mod_time).strftime("%Y-%m-%d")
    print(f"v{__version__} - {mod_time_str}", flush=True)
    print(niceline)


def load_config(config):
    """Import the user-specified config-file to get get OpticalChain(s) and the options-dictionaries.
    The user is reponsible that no harmful code is written in that config-file.
    """
    print("...setting up and importing optical chain(s)...", end="", flush=True)
    # import the required variables from the config-file
    #config = __import__(config_file)

    if hasattr(config, "OpticalChainList"):
        OpticalChainList = config.OpticalChainList
    elif hasattr(config, "OpticalChain"):
        OpticalChainList = config.OpticalChain
    else:
        raise (
            ValueError,
            "Could not import an optical-chain-object or list thereof with the name OpticalChain or OpticalChainList.",
        )

    if hasattr(config, "SourceProperties"):
        SourceProperties = config.SourceProperties
    else:
        print("No SourceProperties-dictionary provided, will use defaults.")
        SourceProperties = {}

    if hasattr(config, "DetectorOptions"):
        DetectorOptions = config.DetectorOptions
    else:
        print("No DetectorOptions-dictionary provided, will use defaults.")
        DetectorOptions = {}

    if hasattr(config, "AnalysisOptions"):
        AnalysisOptions = config.AnalysisOptions
    else:
        print("No AnalysisOptions-dictionary provided, will use defaults.")
        AnalysisOptions = {}

    print(
        "\r\033[K", end="", flush=True
    )  # move to beginning of the line with \r and then delete the whole line with \033[K

    return OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions


def complete_defaults(SourceProperties, DetectorOptions, AnalysisOptions):
    """Load default-options-dictionaries and update them with the values from the user-supplied option-dictionaries.
    Then return those updated dictionaries as the final configuration.
    """

    from ART.DefaultOptions import DefaultSourceProperties, DefaultDetectorOptions, DefaultAnalysisOptions

    DefaultSourceProperties.update(SourceProperties)
    DefaultDetectorOptions.update(DetectorOptions)
    DefaultAnalysisOptions.update(AnalysisOptions)

    return DefaultSourceProperties, DefaultDetectorOptions, DefaultAnalysisOptions


def setup_detector(OpticalChain, DetectorOptions, RayList=None):
    """Set up a detector based on the given parameters.
    Returns a Detector object.
    """

    if DetectorOptions["ManualDetector"]:
        if DetectorOptions["DetectorCentre"] is None:
            raise RuntimeError(
                'For manual detector placement you need to specify "DetectorCentre" in the "DetectorOptions"-dictionary.'
            )
        if DetectorOptions["DetectorNormal"] is None:
            raise RuntimeError(
                'For manual detector placement you need to specify "DetectorNormal" in the "DetectorOptions"-dictionary.'
            )
        Detector = mdet.Detector(
            OpticalChain.optical_elements[DetectorOptions["ReflectionNumber"]].position,
            DetectorOptions["DetectorCentre"],
            DetectorOptions["DetectorNormal"],
        )
    else:
        if DetectorOptions["DistanceDetector"] is None:
            raise RuntimeError(
                'For automatic detector placement you need to specify "DistanceDetector" in the "DetectorOptions"-dictionary.'
            )
        if RayList is None:
            raise RuntimeError(
                'For automatic detector placement you need to add a RayList as an input (selected from the "RayListHistory" by the index DetectorOptions["ReflectionNumber"]).'
            )
        Detector = mdet.Detector(OpticalChain.optical_elements[DetectorOptions["ReflectionNumber"]].position)
        Detector.autoplace(RayList, DetectorOptions["DistanceDetector"])

    return Detector


def optimize_detector(
    RayListAnalysed,
    Detector,
    DetectorOptions,
    verbose=True,
    maxRaystoConsider=1000,
    IntensityWeighted=False,
    Amplitude=None,
    Precision=3,
):
    """Find the optimal detector position, where the optimization-metrix OptFor can be
    "spotsize" (smallest), "duration" (shortest/planest wavefront), or "intensity" (for maximized 1/(spotsize^2*duration)).
    IntensityWeighted specifies whether to weight the optimization by ray-intensity.
    Amplitude specifies the range around the initial detector distance to search for the optimal position.
    Precision specifies the precision of the search algorithm (10^(-precision) fraction of the search Amplitude).
    maxRaystoConsider specifies the maximum number of rays to use for the optimization: if RayListAnalysed contains more rays,
    we'll randomly sample #maxRaystoConsider of them to speed up optimization.

    Returns the the optimal detector distance, optimal detector object (shifted original Detector),
    and the corresponding value of the optimization-metric..
    """
    if len(RayListAnalysed) > maxRaystoConsider:
        RayListForOpt = np.random.choice(RayListAnalysed, maxRaystoConsider, replace=False)
    else:
        RayListForOpt = RayListAnalysed

    OptDetector, OptSizeSpot, OptDuration = mp.FindOptimalDistance(
        Detector, RayListForOpt, DetectorOptions["OptFor"], Amplitude, Precision, IntensityWeighted, verbose
    )
    OptDistance = OptDetector.get_distance()

    if verbose:
        resultstring = f"The optimal detector distance is {OptDistance:.3f} mm, with"
        if IntensityWeighted:
            resultstring += " intensity-weighted"
        if DetectorOptions["OptFor"] in ["intensity", "spotsize"]:
            resultstring += f" spatial std of {OptSizeSpot*1e3:.3g} \u03BCm"
        if DetectorOptions["OptFor"] in ["intensity", "duration"]:
            resultstring += f" temporal std of {OptDuration:.3g} fs."
        print(resultstring, flush=True)

    # Detector = OptDetector.copy_detector()

    return OptDetector, OptSizeSpot, OptDuration


def make_plots(OpticalChain, RayListAnalysed, Detector, SourceProperties, DetectorOptions, AnalysisOptions):
    """Go through the built-in plotting options."""

    if AnalysisOptions["plot_Render"]:
        mplots.RayRenderGraph(
            OpticalChain,
            Detector.get_distance() * 1.2,
            AnalysisOptions["maxRaysToRender"],
            AnalysisOptions["OEPointsToRender"],
            AnalysisOptions["OEPointsScale"],
            draw_mesh = AnalysisOptions["draw_mesh"],
            cycle_ray_colors = AnalysisOptions["cycle_ray_colors"]
        )

    if AnalysisOptions["plot_DelayMirrorProjection"]:
        mplots.MirrorProjection(OpticalChain, DetectorOptions["ReflectionNumber"], Detector, "Delay")

    if AnalysisOptions["plot_IntensityMirrorProjection"]:
        mplots.MirrorProjection(OpticalChain, DetectorOptions["ReflectionNumber"], Detector, "Intensity")

    if AnalysisOptions["plot_IncidenceMirrorProjection"]:
        mplots.MirrorProjection(OpticalChain, DetectorOptions["ReflectionNumber"], Detector, "Incidence")

    if AnalysisOptions["plot_SpotDiagram"]:
        mplots.SpotDiagram(RayListAnalysed, Detector, AnalysisOptions["DrawAiryAndFourier"])

    if AnalysisOptions["plot_DelaySpotDiagram"]:
        mplots.SpotDiagram(RayListAnalysed, Detector, AnalysisOptions["DrawAiryAndFourier"], "Delay")

    if AnalysisOptions["plot_IntensitySpotDiagram"]:
        mplots.SpotDiagram(RayListAnalysed, Detector, AnalysisOptions["DrawAiryAndFourier"], "Intensity")

    if AnalysisOptions["plot_IncidenceSpotDiagram"]:
        mplots.SpotDiagram(RayListAnalysed, Detector, AnalysisOptions["DrawAiryAndFourier"], "Incidence")

    if AnalysisOptions["plot_DelayGraph"]:
        mplots.DelayGraph(
            RayListAnalysed, Detector, SourceProperties["DeltaFT"], AnalysisOptions["DrawAiryAndFourier"], "Delay"
        )

    if AnalysisOptions["plot_IntensityGraph"]:
        mplots.DelayGraph(
            RayListAnalysed, Detector, SourceProperties["DeltaFT"], AnalysisOptions["DrawAiryAndFourier"], "Intensity"
        )

    if AnalysisOptions["plot_IncidenceGraph"]:
        mplots.DelayGraph(
            RayListAnalysed, Detector, SourceProperties["DeltaFT"], AnalysisOptions["DrawAiryAndFourier"], "Incidence"
        )

    if __name__ == "__main__":
        mplots.show()


# %%
def run_ART(OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions, loop=False):
    niceline = "___________________________________________________________________________________________________\n"

    # mp._tic() ################

    """ THE ACTUAL RAYTRACING CALCULATION """
    output_rays = OpticalChain.get_output_rays()
    RayListAnalysed = output_rays[DetectorOptions["ReflectionNumber"]]

    ETransmission = mplots.getETransmission(OpticalChain.source_rays, RayListAnalysed)
    if AnalysisOptions["verbose"]:
        print(niceline[:-1], flush=True)
        if isinstance(OpticalChain.description, str) and len(OpticalChain.description) > 0:
            print("***" + OpticalChain.description + "*** :")
        if OpticalChain.loop_variable_name is not None and OpticalChain.loop_variable_value is not None:
            print(
                "For "
                + OpticalChain.loop_variable_name
                + " = "
                + "{:f}".format(OpticalChain.loop_variable_value)
                + ":\n"
            )
            print("The optical setup has an energy transmission of " + "{:.1f}".format(ETransmission) + "%.\n")

    """ SET UP DETECTOR """
    Detector = setup_detector(OpticalChain, DetectorOptions, RayListAnalysed)

    """ OPTIMIZE DETECTOR POSITION, or just give a RESULT-SUMMARY at the initially specified detector position """
    if DetectorOptions["AutoDetectorDistance"]:
        Detector, SpotSizeSD, DurationSD = optimize_detector(
            RayListAnalysed,
            Detector,
            DetectorOptions,
            AnalysisOptions["verbose"],
            maxRaystoConsider=1000,
            IntensityWeighted=True,
        )
    else:
        SpotSizeSD, DurationSD = mplots.GetResultSummary(Detector, RayListAnalysed, AnalysisOptions["verbose"])

    # mp._toc() ################

    if AnalysisOptions["verbose"]:
        print(niceline)

    """ PRODUCE AND SHOW THE PLOTS AS SELECTED IN AnalysisOptions """
    if not loop or __name__ != "__main__":
        # if any of the keys in AnalysisOptions that start with "plot_" are true; otherwise just do nothing
        plot_keys = [key for key in AnalysisOptions if key.startswith("plot_")]
        if any(AnalysisOptions[key] for key in plot_keys):
            make_plots(OpticalChain, RayListAnalysed, Detector, SourceProperties, DetectorOptions, AnalysisOptions)

    return OpticalChain, Detector, ETransmission, SpotSizeSD, DurationSD


# %%
def main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions, save_file_name=None):
    """COMPLETE CONFIG-DOCTIONARIES WITH DEFAULT OPTIONS"""
    SourceProperties, DetectorOptions, AnalysisOptions = complete_defaults(
        SourceProperties, DetectorOptions, AnalysisOptions
    )

    """ Go through supplied list of OpticalChain(s) and save the results in a dictionary"""
    # initialize a dictionary with the names of the variables to be saved as keys
    keeper_names = ["OpticalChain", "Detector", "ETransmission", "SpotSizeSD", "DurationSD"]
    kept_data = {name: [] for name in keeper_names}

    if isinstance(OpticalChainList, moc.OpticalChain):  # single OpticalChain not in a list
        OpticalChainList = [OpticalChainList]
        loop = False
    elif not isinstance(OpticalChainList, list):
        raise ValueError(
            "The supplied OpticalChain is neither an OpticalChain-object, nor a list of those, as it should be."
        )
    else:
        loop = True

    # do the simulations
    for i, OpticalChain in enumerate(OpticalChainList):
        print("Optical Chain " + str(i) + "/" + str(len(OpticalChainList)) + " ", end="", flush=True)
        OpticalChain, Detector, ETransmission, SpotSizeSD, DurationSD = run_ART(
            OpticalChain, SourceProperties, DetectorOptions, AnalysisOptions, loop
        )
        for name in keeper_names:
            kept_data[name].append(locals()[name])

    # save results if user said so
    if AnalysisOptions["save_results"]:
        print("...saving data...", end="", flush=True)
        mp.save_compressed(kept_data, save_file_name)
        print(
            "\r\033[K", end="", flush=True
        )  # move to beginning of the line with \r and then delete the whole line with \033[K

    return kept_data


# %%
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python ARTmain.py CONFIG_FILE")
        config_file = "CONFIG_test"

        """ LOAD CONFIGURATION """
        OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions = load_config(config_file)

        """ LAUNCH MAIN """
        main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions, save_file_name=config_file)
        # This was simply done to simplify profiling and debugging
        """SHOW THE NAME BANNER"""
        print_banner(1)


        """ LOAD CONFIGURATION """
        OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions = load_config(config_file)

        """ LAUNCH MAIN """
        main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions, save_file_name=config_file)

    else:
        """SHOW THE NAME BANNER"""
        print_banner(1)

        config_file = sys.argv[1]
        filename = os.path.basename(config_file)
        spec = importlib.util.spec_from_file_location(filename, config_file)
        config_module = importlib.util.module_from_spec(spec)
        sys.modules[filename] = config_module
        spec.loader.exec_module(config_module)
        
        """ LOAD CONFIGURATION """
        OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions = load_config(config_module)

        """ LAUNCH MAIN """
        main(OpticalChainList, SourceProperties, DetectorOptions, AnalysisOptions, save_file_name=config_file)
