---
title: "DetectorOptions"
weight: 3
type: docs
---

{{% alert title="Summary" color="primary" %}}
`DetectorOptions` is a python dictionnary
{{% /alert %}}

The detector is a virtual plane of infinite size which plays the role of a CCD. It can be positionned manually by specifying the position of a point on the detector and its normal vector. However, in the case of ART it made sense to make the positionning of the detector part of the main functionality of the code. Thus, the `main()` function accepts the `DetectorOptions` argument. This object allows the code to automatically position the virtual detection plane in the optimal spot. It's a Python dictionnary with the following possible attributes:
 - `"ReflectionNumber"`: This is the number of the optic after the which the detector will be positionned. In most cases that would be the final optic in the setup and we commonly use the Python shorthand for the index of the last element: -1
 - `"ManualDetector"`: Setting this to `True` would disable the automatic positionning of the detector and you would have to also specify `"DetectorCentre"` and `"DetectorNormal"`  to position it.
 - `"AutoDetectorDistance"`: If set to `True`, the algorithm will automatically optimize the detector position based on the `"OptFor"`, `"MaxRaystoConsider"` and `"IntensityWeighted"` options.
 - `"DistanceDetector"`: If `"AutoDetectorDistance"` is `True`, it will be the starting point for the optimisation. Otherwise, this will be the distance of the detector from the last optical element.


```Python
DetectorOptions = {
    'ReflectionNumber' : -1,
    'ManualDetector' : False,
    'DistanceDetector' : 76 ,
    'AutoDetectorDistance' : True,
    'OptFor' : "intensity"
}
```
