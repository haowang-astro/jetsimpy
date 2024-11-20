# define your own radiation model in a pure Python way

import numpy as np
from matplotlib import pyplot as plt
import jetsimpy

# some constants
CSpeed = 29979245800.0
MassP = 1.672622e-24
MassE = 9.109384e-28
SigmaT = 6.6524587e-25
ECharge = 4.803204673e-10

# put the parameters in a dictionary
P = dict(
    Eiso = 1e52,        # Isotropic equivalent energy
    lf = 300,           # Lorentz factor
    theta_c = 0.1,      # half opening angle
    n0 = 1,             # ism number density
    A = 0,              # wind number density amplitude
    eps_e = 0.1,       # epsilon_e
    eps_b = 0.01,      # epsilon_b
    p = 2.17,           # electron power index
    theta_v = 0.4,      # viewing angle (rad)
    d = 474.33,         # distance (Mpc)
    z = 0.1,            # redshift
)

# synchrotron radiation with the deep Newtonian phase described in Sironi & Giannios 2013
def sync_dnp(nu, Psync, blast):
    eps_e = Psync["eps_e"]  # epsilon_e
    eps_b = Psync["eps_b"]  # epsilon_B
    p = Psync["p"]          # electron power index

    n_blast = blast.n_blast # post-shock number density
    t = blast.t             # time since burst
    gamma = blast.gamma     # Lorentz factor of the blast
    e = blast.e_density     # post-shock energy density (rest mass excluded)

    gamma_m = (p - 2.0) / (p - 1.0) * (eps_e * MassP / MassE * (gamma - 1.0))
    f = 1
    if gamma_m <= 1:        # if gamma_m < 1, truncate the electron population
        gamma_m = 1.0
        f = (p - 2.0) / (p - 1.0) * eps_e * MassP / MassE * (gamma - 1.0) / gamma_m
    
    B = np.sqrt(8.0 * np.pi * eps_b * e)
    gamma_c = 6.0 * np.pi * MassE * gamma * CSpeed / SigmaT / B / B / t
    nu_m = 3.0 * ECharge * B * gamma_m * gamma_m / 4.0 / np.pi / CSpeed / MassE
    nu_c = 3.0 * ECharge * B * gamma_c * gamma_c / 4.0 / np.pi / CSpeed / MassE
    e_p = f * np.sqrt(3.0) * ECharge * ECharge * ECharge * B * n_blast / MassE / CSpeed / CSpeed

    if nu_m < nu_c:
        if nu < nu_m:
            emissivity = e_p * np.cbrt(nu / nu_m)
        elif nu < nu_c:
            emissivity = e_p * np.power(nu / nu_m, - (p - 1) / 2.0)
        else:
            emissivity = e_p * np.power(nu_c / nu_m, - (p - 1) / 2.0) * np.power(nu / nu_c, - p / 2)
    else:
        if nu < nu_c:
            emissivity = e_p * np.cbrt(nu / nu_c)
        elif nu < nu_m:
            emissivity = e_p / np.sqrt(nu / nu_c)
        else:
            emissivity = e_p / np.sqrt(nu_m / nu_c) * np.power(nu / nu_m, - p / 2)
    
    isotropic_intensity = emissivity * blast.dR
    return isotropic_intensity

# define the observing time and frequency
tday = np.logspace(-2, 4, 100)
tsecond = tday * 3600 * 24
nu = 3e9

# hydro simulation
jet = jetsimpy.Jet(
    jetsimpy.Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
    0.0,               # wind number density scale
    P["n0"],           # ism number density scale
    tmax=1e11,
    grid=jetsimpy.ForwardJetRes(P["theta_c"], 129)
)

# flux density (default synchrotron model)
flux_density1 = jet.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary
    rtol=1e-3,
    model="sync"
)

# flux density (deep Newtonian phase)
flux_density2 = jet.FluxDensity(
    tsecond,           # [second] observing time span
    nu,                # [Hz]     observing frequency
    P,                 # parameter dictionary
    rtol=1e-3,
    model=sync_dnp
)

# plot the light curves
plt.plot(tday, flux_density1, label="default synchrotron", color="black", linestyle="--")
plt.plot(tday, flux_density2, label="Deep Newtonian Phase", color="black")
plt.xlim(1e-2, 1e4)
plt.ylim(1e-5, 1e1)
plt.xscale("log")
plt.yscale("log")
plt.xlabel("time [day]")
plt.ylabel("flux density [mJy]")
plt.legend()
plt.show()