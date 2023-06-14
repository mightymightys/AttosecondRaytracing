---
title: "CONFIG_2toroidals_f-x-f.py"
weight: 2
---

The config file begins with the necessary imports.

```Python
import numpy as np
import ART.ModuleMirror as mmirror
import ART.ModuleSupport as msupp
import ART.ModuleProcessing as mp
import ART.ModuleMask as mmask
```

```Python
SourceProperties = {
    'Divergence' : 50e-3/2,
    'SourceSize' : 0,
    'Wavelength' : 80e-6,
    'DeltaFT'    : 0.5,
    'NumberRays' : 1000
}
```
