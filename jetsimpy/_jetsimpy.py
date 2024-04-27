import numpy as np
from numpy.typing import NDArray
from typing import Dict
from . import jetsimpy_extension as _extension
from ._grid import Uniform

_C = 29979245800.0
_Mass_P = 1.672622e-24
_mJy = 1e-26
_MPC = 3.09e24
_MAS = 1.0 / 206264806.24709466

class Jet:
    def __init__(             # It is the user's responsibility to make sure input values are valid.
        self,
        theta: NDArray,       # [tabulated data]: polar angles
        energy: NDArray,      # [tabulated data]: Eiso (erg) (rest mass excluded)
        lf: NDArray,          # [tabulated data]: Lorentz factor
        nwind: float,         # [wind density scale]: n = nwind * (r / 1e17)^-2 + nism (cm^-3)
        nism: float,          # [wind density scale]: n = nwind * (r / 1e17)^-2 + nism (cm^-3)
        tmin=10.0,            # [simulation start time]: (s)
        tmax=3.2e9,           # [simulation end time]: (s)
        grid=Uniform(257),    # [cell edge angles]: start with 0 and end with pi.
        tail=True,            # [isotropic tail]: add an extremely low energy low velocity isotropic tail for safty
        spread=True,          # [spreading]: spread or not
        calib_level=1,        # [calibration level]: 0: no calibration. 1: BM all time. 2: smoothly go from BM to ST (dangerous)
        rtol=1e-6,            # [primitive variable solver tolerance]: Don't change it unless you know what is going on.
        cfl=0.9,              # [cfl number]: Don't change it unless you know what is going on.
    ):

        # save
        self.theta_data = theta
        self.energy_data = energy
        self.lf_data = lf
        self.nwind = nwind
        self.nism = nism
        self.tmin = tmin
        self.tmax = tmax
        self.rtol = rtol
        self.cfl = cfl
        self.grid = grid
        self.tail = tail
        self.spread = spread
        self.calib_level = calib_level

        # solve jet
        jet_config = self._configJet()
        self._jet = _extension.Jet(jet_config)
        self._jet.solveJet()

    def _configJet(self):
        # initialize parameter object
        jet_config = _extension.JetConfig()
        jet_config.nwind = self.nwind
        jet_config.nism = self.nism
        jet_config.tmin = self.tmin
        jet_config.tmax = self.tmax
        jet_config.rtol = self.rtol
        jet_config.cfl = self.cfl
        jet_config.spread = self.spread
        jet_config.calib_level = self.calib_level

        # generate grid
        theta_edge = self.grid
        theta = np.array([(theta_edge[i] + theta_edge[i + 1]) / 2 for i in range(len(theta_edge) - 1)])
        jet_config.theta_edge = theta_edge

        # add isotropic tail
        if self.tail:
            self.energy_data[self.energy_data <= np.max(self.energy_data) * 1e-12] = np.max(self.energy_data) * 1e-12
            self.lf_data[self.lf_data <= 1.005] = 1.005

        # interpolate initial condition to grid points
        E0 = np.interp(theta, self.theta_data, self.energy_data / 4.0 / np.pi / _C ** 2.0)
        lf0 = np.interp(theta, self.theta_data, self.lf_data)

        # get mej and initial velocity
        Mej0 = E0 / (lf0 - 1)
        beta0 = np.sqrt(1.0 - 1.0 / lf0 ** 2)
        
        # analytically expand (coast) the blastwave from t=0 to t=tmin
        R0 = beta0 * _C * self.tmin
        Msw0 = self.nwind * _Mass_P * R0 / 1e17 * 1e51 + self.nism * _Mass_P * R0 * R0 * R0 / 3.0
        Eb0 = E0 + Mej0 + Msw0
        
        # config jet initial condition
        jet_config.Eb = Eb0
        jet_config.Ht = np.zeros_like(theta)
        jet_config.Msw = Msw0
        jet_config.Mej = Mej0
        jet_config.R = R0
        
        return jet_config
    
    # ---------- PDE original data ---------- #
    @property
    def t_pde(self) -> NDArray:
        result = np.array(self._jet.getT())
        return result
    
    # original y: Msw, Mej, beta_gamma_sq, beta_th, R (shape = [5, ntheta, nt])
    @property
    def y_pde(self) -> NDArray:
        result = np.array(self._jet.getY())
        return result
    
    @property
    def theta_pde(self) -> NDArray:
        result = np.array(self._jet.getTheta())
        return result

    # ---------- PDE data interpolation ---------- #
    # rest mass excluded energy (erg / sr)
    def dE0_dOmega(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateE0(t, theta)
    
    # swetp-up mass (g / sr)
    def dMsw_dOmega(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateMsw(t, theta)
    
    # ejecta mass (g / sr)
    def dMej_dOmega(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateMej(t, theta)

    # four velocity
    def beta_gamma(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateBetaGamma(t, theta)
    
    # polar velocity
    def beta_theta(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateBetaTh(t, theta)
    
    # radius (cm)
    def R(self, t: NDArray, theta: NDArray) -> NDArray:
        return self._jet.interpolateR(t, theta)

    # ---------- Radiation Related ---------- #

    # specific intensity at jet sphreical coordinate [cgs]
    def Intensity(self, t: NDArray, nu: NDArray, theta: NDArray, phi: NDArray, para: Dict[str, float], model="sync") -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        try:
            I = self._jet.calculateIntensity(t, nu, theta, phi)
        except Exception as e:
            raise e
        
        return I
    
    # flux density [mJy]
    def FluxDensity(self, t: NDArray, nu: NDArray, para: Dict[str, float], model="sync", rtol=1e-3) -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        try:
            L = self._jet.calculateLuminosity(t, nu, rtol)
        except Exception as e:
            raise e
        
        return L * (1 + para["z"]) / 4 / np.pi / (para["d"] * _MPC) ** 2 / _mJy
    
    # apparent superluminal motion [mas]
    def Offset(self, t: NDArray, nu: NDArray, para: Dict[str, float], model="sync", rtol=1e-3) -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        # config offset model
        self._jet.configAvgModel("offset")

        # calculate weighted average
        try:
            integral = self._jet.calculateAvgModel(t, nu, rtol)
            luminosity = self._jet.calculateLuminosity(t, nu, rtol)
        except Exception as e:
            raise e
        weighted = integral / luminosity
        
        return weighted / para["d"] / _MPC / (1.0 + para["z"]) / (1.0 + para["z"]) / _MAS
    
    # size along the jet axis [mas]
    def SizeX(self, t: NDArray, nu: NDArray, para: Dict[str, float], model="sync", rtol=1e-3) -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        # ---------- calculate offset ---------- #

        # config offset model
        self._jet.configAvgModel("offset")

        # calculate weighted average
        try:
            integral_offset = self._jet.calculateAvgModel(t, nu, rtol)
            luminosity = self._jet.calculateLuminosity(t, nu, rtol)
        except Exception as e:
            raise e
        xc = integral_offset / luminosity

        # ---------- calculate xscale ---------- #

        # config x scale model
        self._jet.configAvgModel("sigma_x")

        # calculate weighted average
        try:
            integral_x_sq = self._jet.calculateAvgModel(t, nu, rtol)
            luminosity = self._jet.calculateLuminosity(t, nu, rtol)
        except Exception as e:
            raise e
        x_sq = integral_x_sq / luminosity

        # calculate xscale. The following notes are for myself in case I forget what is going on.
        # This makes sense because of the following reason:
        # First, ∫x dL = xc * ∫dL based on xc defination.
        # Therefore, 
        # sigma^2_x = ∫(x - xc)^2 dL / ∫dL 
        #           = ∫(x^2 - 2*x*xc + xc^2)dL / ∫dL 
        #           = (∫x^2 dL - 2*xc*∫x dL + xc^2*∫dL) / ∫dL
        #           = (∫x^2 dL - 2*xc^2*∫dL + xc^2*∫dL) / ∫dL
        #           = (∫x^2 dL / ∫dL) - xc^2
        sigma_x = np.sqrt(x_sq - xc * xc)
        
        return sigma_x / para["d"] / _MPC / (1.0 + para["z"]) / (1.0 + para["z"]) / _MAS
    
    # size perpendicular to the jet axis [mas]
    def SizeY(self, t: NDArray, nu: NDArray, para: Dict[str, float], model="sync", rtol=1e-3) -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        # config y scale model
        self._jet.configAvgModel("sigma_y")

        # calculate weighted average
        try:
            integral_y_sq = self._jet.calculateAvgModel(t, nu, rtol)
            luminosity = self._jet.calculateLuminosity(t, nu, rtol)
        except Exception as e:
            raise e
        
        # calculate y scale
        sigma_y = np.sqrt(integral_y_sq / luminosity)
        
        return sigma_y / para["d"] / _MPC / (1.0 + para["z"]) / (1.0 + para["z"]) / _MAS

    # [cgs] specific intensity at LOS frame coordinate (x_tilde, y_tilde). This method is intended for sky map.
    def IntensityOfPixel(self, t: NDArray, nu: NDArray, x_tilde: NDArray, y_tilde: NDArray, para: Dict[str, float], model="sync") -> NDArray:
        # config parameters
        self._jet.configParameters(para)

        # config emissivity model
        self._jet.configEmissivity(model)

        return self._jet.IntensityOfPixel(t, nu, x_tilde, y_tilde)