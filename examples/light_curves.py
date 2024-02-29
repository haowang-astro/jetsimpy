# light curves of a gaussian jet

import numpy as np
from matplotlib import pyplot as plt
import jetsimpy

P = dict(
    Eiso = 1e52,        # (Jet) Isotropic equivalent energy
    lf = 300,           # (Jet) Lorentz factor
    theta_c = 0.1,      # (Jet) half opening angle
    n0 = 1,             # (ISM) constant number density
    k = 2.0,            # (ISM) wind power index
    A = 0,              # (ISM) wind amplitude
    eps_e = 0.1,        # (Radiation) epsilon_e
    eps_b = 0.01,       # (Radiation) epsilon_b
    p = 2.17,           # (Radiation) electron power index
    theta_v = 0.4,      # (Radiation) viewing angle
    d = 474.33,         # (radiation) distance (Mpc)
    z = 0.1,            # (radiation) redshift
    b = 0,              # (radiation) magnetic field anisotropy
)

# tabulated energy/LF structure
theta = np.linspace(0, np.pi, 1000)
Eiso = P["Eiso"] * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2)
lf = (P["lf"] - 1) * np.exp(- 0.5 * (theta / P["theta_c"]) ** 2) + 1

# solve the evolution of the jet
# jet without spreading
jet1 = jetsimpy.Afterglow(
    theta,           # array of theta
    Eiso,            # array of isotropic equivalent energy
    lf,              # array of initial lorentz factor
    P["A"],          # scale of wind density
    P["n0"],         # constant number density
    spread=False,    # (default = True) with/without spreading effect 
    coast=True,      # (default = True) with/without coasting. If this is "False", the initial lorentz factor data will be omitted.
)

# jet with spreading
jet2 = jetsimpy.Afterglow(
    theta,           # array of theta
    Eiso,            # array of isotropic equivalent energy
    lf,              # array of initial lorentz factor
    P["A"],          # scale of wind density
    P["n0"],         # constant number density
    spread=True,     # (default = True) with/without spreading effect 
    coast=True,      # (default = True) with/without coasting. If this is "False", the initial lorentz factor data will be omitted.
)

# define the observing time and frequency
tday = np.logspace(-2, 3, 100)
tsecond = tday * 3600 * 24
nu = 3e9

# calculate the afterglow flux density (unit: mJy)
flux1 = jet1.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
    rtol=1e-2,         # (default=1e-2) integration error tolerance
    model="sync",      # default radiation model
)

flux2 = jet2.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
    rtol=1e-2,         # (default=1e-2) integration error tolerance
    model="sync",      # default radiation model
)

# plot the light curves
plt.plot(tday, flux1, label="without spreading", color="black", linestyle="--")
plt.plot(tday, flux2, label="with spreading", color="black")
plt.xlim(1e-2, 1e3)
plt.ylim(1e-5, 1e1)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("time [day]")
plt.ylabel("flux density [mJy]")
plt.legend()
plt.show()