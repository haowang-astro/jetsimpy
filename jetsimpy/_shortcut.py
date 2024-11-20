from ._jetsimpy import Jet
from ._grid import *
from ._jet_type import *

def FluxDensity_tophat(t, nu, P, tmin=10.0, tmax=1e10, spread=True, cal_level=1, cfl=0.9, model="sync", rtol=1e-3, max_iter=100, force_return=True):
    # simulation
    jet = Jet(
        TopHat(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
        P["A"],                        # wind number density scale
        P["n0"],                       # ism number density scale
        tmin=tmin,                     # [simulation start time]: (s)
        tmax=tmax,                    # [simulation end time]: (s)
        grid=ForwardJetRes(P["theta_c"], 129),    # [cell edge angles]: must start with 0 and end with pi.
        tail=True,                     # [isotropic tail]: add an extremely low energy low velocity isotropic tail for safty
        spread=spread,                   # w/wo spreading effect 
        cal_level=cal_level,                   # [calibration level]: 0: no calibration. 1: BM all time. 2: smoothly go from BM to ST (dangerous)
        rtol=1e-6,                     # [primitive variable solver tolerance]: Don't change it unless you know what is going on.
        cfl=cfl,                       # [cfl number]: Don't change it unless you know what is going on.
    )

    # flux density
    flux = jet.FluxDensity(
        t,           # [second] observing time span
        nu,                # [Hz]     observing frequency
        P,                 # parameter dictionary
        model=model,      # emissivity model
        rtol=rtol,         # integration tolerance
        max_iter=max_iter,
        force_return=force_return
    )

    return flux

def FluxDensity_gaussian(t, nu, P, tmin=10.0, tmax=1e10, spread=True, cal_level=1, cfl=0.9, model="sync", rtol=1e-3, max_iter=100, force_return=True):
    # simulation
    jet = Jet(
        Gaussian(P["theta_c"], P["Eiso"], lf0=P["lf"]),    # jet profile
        P["A"],                        # wind number density scale
        P["n0"],                       # ism number density scale
        tmin=tmin,                     # [simulation start time]: (s)
        tmax=tmax,                    # [simulation end time]: (s)
        grid=ForwardJetRes(P["theta_c"], 129),    # [cell edge angles]: must start with 0 and end with pi.
        tail=True,                     # [isotropic tail]: add an extremely low energy low velocity isotropic tail for safty
        spread=spread,                   # w/wo spreading effect 
        cal_level=cal_level,                   # [calibration level]: 0: no calibration. 1: BM all time. 2: smoothly go from BM to ST (dangerous)
        rtol=1e-6,                     # [primitive variable solver tolerance]: Don't change it unless you know what is going on.
        cfl=cfl,                       # [cfl number]: Don't change it unless you know what is going on.
    )

    # flux density
    flux = jet.FluxDensity(
        t,           # [second] observing time span
        nu,                # [Hz]     observing frequency
        P,                 # parameter dictionary
        model=model,      # emissivity model
        rtol=rtol,         # integration tolerance
        max_iter=max_iter,
        force_return=force_return
    )

    return flux

def FluxDensity_powerlaw(t, nu, P, tmin=10.0, tmax=1e10, spread=True, cal_level=1, cfl=0.9, model="sync", rtol=1e-3, max_iter=100, force_return=True):
    # simulation
    jet = Jet(
        PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"], s=P["s"]),    # jet profile
        P["A"],                        # wind number density scale
        P["n0"],                       # ism number density scale
        tmin=tmin,                     # [simulation start time]: (s)
        tmax=tmax,                    # [simulation end time]: (s)
        grid=ForwardJetRes(P["theta_c"], 129),    # [cell edge angles]: must start with 0 and end with pi.
        tail=True,                     # [isotropic tail]: add an extremely low energy low velocity isotropic tail for safty
        spread=spread,                   # w/wo spreading effect 
        cal_level=cal_level,                   # [calibration level]: 0: no calibration. 1: BM all time. 2: smoothly go from BM to ST (dangerous)
        rtol=1e-6,                     # [primitive variable solver tolerance]: Don't change it unless you know what is going on.
        cfl=cfl,                       # [cfl number]: Don't change it unless you know what is going on.
    )

    # flux density
    flux = jet.FluxDensity(
        t,           # [second] observing time span
        nu,                # [Hz]     observing frequency
        P,                 # parameter dictionary
        model=model,      # emissivity model
        rtol=rtol,         # integration tolerance
        max_iter=max_iter,
        force_return=force_return
    )

    return flux
