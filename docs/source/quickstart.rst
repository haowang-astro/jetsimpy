Quickstart
==========

Analytic jet profiles
---------------------
The program supports an easy way to calculate the light curve or spectrum if the jet profile is Top-hat, Gaussian, or Power-law. 

For Top-hat jets::

    # put the parameters in a dictionary
    P = dict(
        Eiso = 1e52,        # Core isotropic equivalent energy
        lf = 300,           # Core Lorentz factor
        theta_c = 0.1,      # half opening angle
        n0 = 1,             # ism number density
        A = 0,              # wind number density amplitude
        eps_e = 0.1,        # epsilon_e
        eps_b = 0.01,       # epsilon_b
        p = 2.17,           # electron power index
        theta_v = 0.0,      # viewing angle (rad)
        d = 474.33,         # distance (Mpc)
        z = 0.1,            # redshift
    )

    # define the observing time and frequency
    tday = np.logspace(-3, 3, 100)
    tsecond = tday * 3600 * 24
    nu = 1e18

    # flux density
    fd_tophat = jetsimpy.FluxDensity_tophat(tsecond, nu, P)

For Gaussian jets::

    # put the parameters in a dictionary
    P = dict(
        Eiso = 1e52,        # Core isotropic equivalent energy
        lf = 300,           # Core Lorentz factor
        theta_c = 0.1,      # half opening angle
        n0 = 1,             # ism number density
        A = 0,              # wind number density amplitude
        eps_e = 0.1,        # epsilon_e
        eps_b = 0.01,       # epsilon_b
        p = 2.17,           # electron power index
        theta_v = 0.0,      # viewing angle (rad)
        d = 474.33,         # distance (Mpc)
        z = 0.1,            # redshift
    )

    # define the observing time and frequency
    tday = np.logspace(-3, 3, 100)
    tsecond = tday * 3600 * 24
    nu = 1e18

    # flux density
    fd_gaussian = jetsimpy.FluxDensity_gaussian(tsecond, nu, P)

For power-law jets::

    # put the parameters in a dictionary
    P = dict(
        Eiso = 1e52,        # Core isotropic equivalent energy
        lf = 300,           # Core Lorentz factor
        theta_c = 0.1,      # half opening angle
        n0 = 1,             # ism number density
        A = 0,              # wind number density amplitude
        eps_e = 0.1,        # epsilon_e
        eps_b = 0.01,       # epsilon_b
        p = 2.17,           # electron power index
        theta_v = 0.0,      # viewing angle (rad)
        d = 474.33,         # distance (Mpc)
        z = 0.1,            # redshift
        s = 6,              # power-law jet slope
    )

    # define the observing time and frequency
    tday = np.logspace(-3, 3, 100)
    tsecond = tday * 3600 * 24
    nu = 1e18

    # flux density
    fd_powerlaw = jetsimpy.FluxDensity_powerlaw(tsecond, nu, P)

Tabulated jet profiles
----------------------
For tabulated jet profiles, you will need two steps to calculate the GRB observables. The first step is to simulate the hydrodynamic evolution of a jet. The second step is to calculate the afterglow observables on top of the simulation results.

Below is an example::

  import numpy as np
  import jetsimpy

  # generate some tabulated jet profiles
  theta = np.linspace(0, np.pi, 10000)                      # polar angles
  Eiso = 1e53 * np.exp(- 0.5 * (theta / 0.1) ** 2)          # isotropic equivalent energy
  lf = (1000 - 1) * np.exp(- 0.5 * (theta / 0.1) ** 2) + 1  # initial Lorentz factor

  # parameter dictionary
  P = dict(
      eps_e = 0.1,        # epsilon_e
      eps_b = 0.01,       # epsilon_b
      p = 2.17,           # electron power index
      theta_v = 0.5,      # viewing angle (rad)
      d = 474.33,         # luminosity distance (Mpc)
      z = 0.1,            # redshift
  )

  # (step 1) hydro simulation
  jet = jetsimpy.Jet(
      (theta, Eiso, lf), # specify the jet profile here
      0.0,               # wind number density scale
      1,                 # ism number density scale
      grid=jetsimpy.ForwardJetRes(0.1, 129),    # simulation resolution
  )
  
  # (step 2) flux density [mJy]
  tday = np.logspace(-2, 3, 100)  # observing time in day
  tsecond = tday * 3600 * 24      # observing time in second
  nu = 1e15                       # observing frequency

  flux_density = jet.FluxDensity(
      tsecond,           # [second] observing time span
      nu,                # [Hz]     observing frequency
      P,                 # parameter dictionary
  )

Tips on the resolution
----------------------
The resolution setup is specified by the ``grid`` keyword which takes an array ranging from 0 to :math:`\pi`.

For any jet profile, the resolution is always **extremely important**. On one hand, if the specified resolution is too sparse to reveal the details of the jet profile evolution, the simulation will not converge. On the other hand, if the resolution is too dense, the simulation will be very time costly. 
An optimal setup is to specify high resolution around the jet core and relatively low resolution outside the core. 

For a jet profile where a half opening angle is well-defined, `jetsimpy` provides the following resolution setups which have been widely tested for both reliability and speed::

    jetsimpy.ForwardJetRes(theta_c, n)         # for forward jet
    jetsimpy.CounterJetRes(theta_c, n)         # for counter jet
    jetsimpy.ForwardCounterJetRes(theta_c, n)  # for both forward and counter jet

where ``n`` is the number of grid points. Details of the resolution can be found in the next section.

For arbitrary jet profiles, unfortunately, there is no universal solution. The default resolution is uniform, but remember this setup is in general **NOT** optimal. It is strongly recomended that users specify their own resolution.