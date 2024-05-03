jetsimpy
===================================
*jetsimpy* is a code for modeling Gamma-ray burst afterglows with arbitrary angular energy and Lorentz factor profile. It includes the hydrodynamic simulation of the relativistic jet, and the synchrotron radiation on top of it.

Installation
===================================
The code is not published to PyPI, so please install it from source::

$ pip install . && python setup.py clean

Quickstart
===================================
This is a simple example of a Gaussian jet::

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
      jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),
      0.0,               # wind number density scale
      P["n0"],           # ism number density scale
  )
  
  # flux density [mJy]
  flux_density = jet.FluxDensity(
      tsecond,           # [second] observing time span
      nu,                # [Hz]     observing frequency
      P,                 # parameter dictionary
  )

More examples can be found in the example folder.