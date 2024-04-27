# jetsimpy
Hydrodynamic simulations of relativistic blastwave with tabulated angular energy and Lorentz factor profiles. Efficient Gamma-Ray Burst afterglow modeling.

Code paper: [Wang et al. 2024](https://arxiv.org/abs/2402.19359).

## Important changes (4. 26. 2024)
I mostly rewrite the package for better code readability, better user interface, better code structure, and (perhaps) better performance. Many variables and functions are renamed, so make sure you update the code to the latest version.

I will try my best to keep the interface unchanged from now on.

## Installation 
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install jetsimpy.
```bash
pip install .
```

Clean up the installation directory. This help avoid conflicts for future installations.
```bash
python setup.py clean
```

You may install the package and clean up the directory with one line:
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
    *jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
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

More examples are available in the example folder. Detailed documentation is being developed.