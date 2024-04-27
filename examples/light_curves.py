# light curves of a gaussian jet

import numpy as np
from matplotlib import pyplot as plt
import jetsimpy

# put the parameters in a dictionary
P = dict(
    Eiso = 1e52,        # (Jet) Isotropic equivalent energy
    lf = 300,           # (Jet) Lorentz factor
    theta_c = 0.1,      # (Jet) half opening angle
    n0 = 1,             # (ISM) constant number density
    A = 0,              # (ISM) wind amplitude
    eps_e = 0.1,        # (Radiation) epsilon_e
    eps_b = 0.01,       # (Radiation) epsilon_b
    p = 2.17,           # (Radiation) electron power index
    theta_v = 0.4,      # (Radiation) viewing angle
    d = 474.33,         # (radiation) distance (Mpc)
    z = 0.1,            # (radiation) redshift
)

# ---------- (step 1) hydro simulation of the jet ---------- #

# jet without spreading
jet1 = jetsimpy.Jet(
    *jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    P["A"],          # scale of wind density
    P["n0"],         # constant number density
    spread=False,    # w/wo spreading effect 
)

# jet with spreading
jet2 = jetsimpy.Jet(
    *jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    P["A"],          # scale of wind density
    P["n0"],         # constant number density
    spread=True,     # w/wo spreading effect 
)

# ---------- (step 2) calculate flux density ---------- #

# define the observing time and frequency
tday = np.logspace(-2, 3, 100)
tsecond = tday * 3600 * 24
nu = 3e9

# calculate the afterglow flux density (unit: mJy)
flux1 = jet1.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
)

flux2 = jet2.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
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