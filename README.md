# jetsimpy
Hydrodynamic simulations of relativistic blastwave with tabulated angular energy and Lorentz factor profiles. Efficient Gamma-Ray Burst afterglow modeling.

Code paper: [Wang et al. 2024](https://arxiv.org/abs/2402.19359).

## Known issues
Sometimes the code fails with an error information "`f(a) and f(b) must have different signs!`" or "`Solver doesn't converge after 100 iterations!`". This means the hydro simulation runs into numerical instability. For now, a temporary solution is to reduce the `cfl` keyword (e.g., 0.5 or smaller). This will slow down the code by a corresponding ratio but make it more stable. This issue will be fixed in the next update.

## Features

### These features are currently supported:
* Tabulated angular energy profile
* Tabulated angular Lorentz factor profile
* ISM / wind / mixed external density profile: $n=n_{\rm ism}+n_{\rm wind}(r/10^{17}{\rm cm})^{-2}$
* Synthetic afterglow light curves
* Apparent superluminal motion
* Sky map and Gaussian equivalent image size

Additionally, you can add your own emissivity model by defining a lambda function in a [c++ source file](jetsimpy/src/Afterglow/models.cpp). This might be helpful if you have a more complicated model such as Synchrotron self-absorption or other cool stuffs. After adding your own model, just install the package as usual and refer to this model with its model name from Python side.

### These features are not supported yet:
* Reverse shock
* Energy injection

## Installation 
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install jetsimpy.
```bash
pip install .
```

Clean up the installation directory. This help avoid conflicts for future installations.
```bash
python setup.py clean
```

You may install the package and clean up the directory by a single line:
```bash
pip install . && python setup.py clean
```

## Quickstart
```Python
import numpy as np
import jetsimpy

# parameter dictionary
P = dict(
    Eiso = 1e53,        # Isotropic equivalent energy (erg)
    lf = 600,           # initial Lorentz factor
    theta_c = 0.1,      # half opening angle
    n0 = 1,             # ism number density
    eps_e = 0.1,        # epsilon_e
    eps_b = 0.01,       # epsilon_b
    p = 2.17,           # electron power index
    theta_v = 0.4,      # viewing angle
    d = 474.33,         # luminosity distance (Mpc)
    z = 0.1,            # redshift
)

# time and frequency
tday = np.logspace(-2, 3, 100)
tsecond = tday * 3600 * 24
nu = 1e15

# hydro simulation
jet = jetsimpy.Jet(
    jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    0.0,               # wind number density scale
    P["n0"],           # ism number density scale
)

# flux density [mJy]
flux_density = jet.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary
)
```

More examples are available in the example folder. 