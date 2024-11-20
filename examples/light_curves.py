# light curves of a gaussian jet

import numpy as np
from matplotlib import pyplot as plt
import jetsimpy

# put the parameters in a dictionary
P = dict(
    Eiso = 1e52,        # Isotropic equivalent energy
    lf = 300,           # Lorentz factor
    theta_c = 0.1,      # half opening angle
    n0 = 1,             # ism number density
    A = 0,              # wind number density amplitude
    eps_e = 0.1,        # epsilon_e
    eps_b = 0.01,       # epsilon_b
    p = 2.17,           # electron power index
    theta_v = 0.4,      # viewing angle (rad)
    d = 474.33,         # distance (Mpc)
    z = 0.1,            # redshift
)

# ---------- (step 1) hydro simulation of the jet ---------- #

# jet without spreading
jet1 = jetsimpy.Jet(
    jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    P["A"],                        # wind number density scale
    P["n0"],                       # ism number density scale
    spread=False,                  # w/wo spreading effect 
    grid=jetsimpy.NorthPole(P["theta_c"], 129)  # resolution
)

# jet with spreading (show full argument and keyword list)
theta, Eiso, lf = jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"])
jet2 = jetsimpy.Jet(
    (theta, Eiso, lf),             # [tuple of tabulated data]: (polar angles, rest mass excluded energy, Lorentz factor)
    P["A"],                        # [wind density scale]: n = nwind * (r / 1e17)^-2 + nism (cm^-3)
    P["n0"],                       # [ism density scale]: n = nwind * (r / 1e17)^-2 + nism (cm^-3)
    tmin=10.0,                     # [simulation start time]: (s)
    tmax=3.2e9,                    # [simulation end time]: (s)
    grid=jetsimpy.NorthPole(P["theta_c"], 129),    # [cell edge angles]: must start with 0 and end with pi.
    tail=True,                     # [isotropic tail]: add an extremely low energy low velocity isotropic tail for safty
    spread=True,                   # w/wo spreading effect 
    cal_level=1,                   # [calibration level]: 0: no calibration. 1: BM all time. 2: smoothly go from BM to ST (dangerous)
    rtol=1e-6,                     # [primitive variable solver tolerance]: Don't change it unless you know what is going on.
    cfl=0.9,                       # [cfl number]: Don't change it unless you know what is going on.
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
    P,                 # parameter dictionary
)

# show full keyword list
flux2 = jet2.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary
    model="sync",      # emissivity model
    rtol=1e-3,         # integration tolerance
    max_iter=100,
    force_return=True
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