import numpy as np
from matplotlib import pyplot as plt
import jetsimpy

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
    s = 6,              # power-law jet slope (required for power-law jet)
)

# define the observing time and frequency
tday = np.logspace(-3, 3, 100)
tsecond = tday * 3600 * 24
nu = 1e18

# flux density
fd_tophat = jetsimpy.FluxDensity_tophat(tsecond, nu, P)
fd_gaussian = jetsimpy.FluxDensity_gaussian(tsecond, nu, P)
fd_powerlaw = jetsimpy.FluxDensity_powerlaw(tsecond, nu, P)

plt.plot(tday, fd_tophat, label="Top-Hat")
plt.plot(tday, fd_gaussian, label="Gaussian")
plt.plot(tday, fd_powerlaw, label="Power-law")

plt.xscale("log")
plt.yscale("log")
plt.xlabel("t [day]")
plt.ylabel(r"$F_\nu$ [mJy]")
plt.legend()
plt.show()