# jetsimpy
Hydrodynamic simulations of relativistic blastwave with tabulated angular energy and Lorentz factor profiles. Efficient Gamma-Ray Burst afterglow modeling.

Code paper: [Wang et al. 2024](https://arxiv.org/abs/2402.19359).

## Updates
Big updates:
* One can now control the maximum iteration number for flux density integration and enforce the program to return a meaningful value even if the integration has not achieved the desired accuracy. This way one can significantly reduce the amount of runtime errors which frequently interrupt the program.
* To define your own radiation model, you now need to define the "solid angle integrated specific intensity" instead of "emissivity" in the previous version. This way it makes more sense to introduce optical depth in the model. The relation is $I_\nu=\frac{j_\nu}{\alpha_\nu}(1-e^{-\alpha_\nu\Delta R})$. For optically thin case ($\alpha_\nu\sim 0$) this relation reduces to $I_\nu=j_\nu\Delta R$.

Small updates:
* For analytic jet profiles (top-hat, Gaussian, power-law), the process to generate synthetic light curves or spectrum is significantly simplified. See "examples/quick_start.py" for more details.
* Some updates in the documentation. The importance of "resolution" is now highlighted.
* The built-in resolution setups are renamed to increase readability. Previous names are still supported for now.

## Documentation
A brief documentation is now available at: [https://jetsimpy.readthedocs.io/](https://jetsimpy.readthedocs.io/).

## Features

### These features are currently supported:
* Tabulated angular energy profile
* Tabulated angular Lorentz factor profile
* ISM / wind / mixed external density profile: $n=n_{\rm ism}+n_{\rm wind}(r/10^{17}{\rm cm})^{-2}$
* Synthetic afterglow light curves
* Apparent superluminal motion
* Sky map and Gaussian equivalent image size

Additionally, you can add your own radiation model by defining a lambda function in a [c++ source file](jetsimpy/src/Afterglow/models.cpp). This might be helpful if you have a more complicated model such as Synchrotron self-absorption or other cool stuffs. After adding your own model, just install the package as usual and refer to this model with its model name from Python side.

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
    A = 0,              # wind number density amplitude
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

# flux density
fd_gaussian = jetsimpy.FluxDensity_gaussian(tsecond, nu, P)
```

More examples are available in the example folder. 