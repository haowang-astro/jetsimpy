# Apparent superluminal motion of a gaussian jet

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

# ---------- (step 2) calculate centroid motion ---------- #

# define the observing time and frequency
tday = np.logspace(-2, 3, 100)
tsecond = tday * 3600 * 24
nu = 3e9

# calculate the afterglow centroid motion (unit: mas)
offset1 = jet1.Offset(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
)

offset2 = jet2.Offset(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary for radiation
)

# calculate the speed of light motion
MPC = 3.086e24
MAS = 4.848e-9
offset_c = 3e10 * tsecond / P["d"] / MPC / (1.0 + P["z"]) / (1.0 + P["z"]) / MAS

# plot the light curves
plt.plot(tday, offset1, label="without spreading", color="black", linestyle="--")
plt.plot(tday, offset2, label="with spreading", color="black")
plt.plot(tday, offset_c, label="speed of light", color="black", linestyle=":")
plt.xlim(0, 1e3)
plt.ylim(0, 0.3)
plt.xlabel("time [day]")
plt.ylabel("Centroid Offset [mas]")
plt.legend()
plt.show()