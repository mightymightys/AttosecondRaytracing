# -*- coding: utf-8 -*-
"""
Created in Jan 2023

This is the entry point for the Attosecond Ray Tracing (ART) software.

@author: Stefan Haessler + André Kalouguine
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
from importlib.metadata import version


def print_banner(i):
    """
    Show the important name banner and the version-number and date.
    """

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
    print("ART version " + version("AttosecondRayTracing"), flush=True)
    print("ARTcore version " + version("AttosecondRayTracing_core"), flush=True)
    print(niceline)


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


# %% Entry point for the script
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
