import numpy as np
from . import jetsimpy_extension
from .autorefine import AutoGrid

C = 29979245800.0
Mass_P = 1.672622e-24
mJy = 1e-26

def LinearInterpolate(xdata, ydata):
    def func(x):
        return np.interp(x, xdata, ydata)
    return func

class Afterglow:
    def __init__(self,
                 theta, Eiso, lf, nwind, nism, 
                 tmin=10, tmax=3.2e9, rtol=1e-6, cfl=0.9, grid=AutoGrid(128), tail=True, spread=True, coast=True, calibration=True):

        self.theta_data = theta
        self.energy_data = Eiso
        self.lf_data = lf
        self.nwind = nwind
        self.nism = nism
        self.rtol = rtol
        self.tmin = tmin
        self.tmax = tmax
        self.cfl = cfl
        self.spread = spread
        self.afterglow = jetsimpy_extension.afterglow()
        self.tail = tail
        self.coast = coast
        self.calibration = calibration

        # preprocess
        E0 = Eiso / C ** 2
        lf0 = lf / 1.0

        # add isotropic tail
        if self.tail:
            E0[E0 <= np.max(E0) * 1e-12] = np.max(E0) * 1e-12
            lf0[lf0 <= 1.005] = 1.005
        
        # interpolate energy and LF
        self.energy = LinearInterpolate(theta, E0 / 4.0 / np.pi)
        self.lf = LinearInterpolate(theta, lf0)

        # generate grid
        if type(grid) == AutoGrid:
            if self.coast:
                self.theta_edge = grid.grid2(self.theta_data, self.energy_data, self.lf_data)
            else:
                self.theta_edge = grid.grid1(self.theta_data, self.energy_data)
        elif type(grid) == np.ndarray:
            self.theta_edge = grid
        else:
            raise RuntimeError("grid type unsupported.")
        self.theta = np.array([(self.theta_edge[i] + self.theta_edge[i + 1]) / 2 for i in range(len(self.theta_edge) - 1)])

        # calibration
        if self.calibration:
            self.afterglow.calibrate_jet()

        # solve jet
        self.solve_jet()
        self.interpolator = self.afterglow.get_interpolator()

    # ---------- Solve the jet dynamics from c++ extension ---------- #
    def solve_jet(self):
        # coasting ?
        if (not self.coast):
            # setup initial condition
            E0 = self.energy(self.theta)
            R0 = C * self.tmin * np.ones_like(self.theta)
            Msw0 = self.nwind * Mass_P * R0 / 1e17 * 1e51 + self.nism * Mass_P * R0 * R0 * R0 / 3.0
            Mej0 = np.zeros_like(self.theta)
            E0 += Msw0
            P0 = np.zeros_like(self.theta)
            y0 = np.array([E0, P0, Msw0, Mej0, R0])

            # config jet
            self.afterglow.config_jet(self.theta_edge, y0, self.spread, self.tmin, self.tmax, self.rtol, self.cfl, self.nwind, self.nism)

            # solve jet
            self.afterglow.solve_jet()
        else:
            # setup initial condition
            E0 = self.energy(self.theta)
            gamma0 = self.lf(self.theta)
            beta0 = np.sqrt(1 - 1 / gamma0 / gamma0)
            betaf0 = 4.0 * beta0 * gamma0 * gamma0 / (4.0 * gamma0 * gamma0 - 1.0)
            R0 = betaf0 * C * self.tmin * np.ones_like(self.theta)
            Mej0 = E0 / (gamma0 - 1)
            Msw0 = self.nwind * Mass_P * R0 / 1e17 * 1e51 + self.nism * Mass_P * R0 * R0 * R0 / 3.0
            E0 += (Msw0 + Mej0)
            P0 = np.zeros_like(self.theta)
            y0 = np.array([E0, P0, Msw0, Mej0, R0])

            # config jet
            self.afterglow.config_jet(self.theta_edge, y0, self.spread, self.tmin, self.tmax, self.rtol, self.cfl, self.nwind, self.nism)

            # solve jet
            self.afterglow.solve_jet()
    
    # ---------- PDE original data ---------- #
    @property
    def t_pde(self):
        return np.array(self.afterglow.get_t_pde(), copy=False)
    
    @property # original y: E, HtR, Msw, Mej, R (shape = [5, ntime, ntheta])
    def y_pde(self):
        return np.array(self.afterglow.get_y_pde(), copy=False).swapaxes(0, 1)

    # ---------- PDE data interpolation ---------- #
    def Eb0(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.Eb0(t, theta)).T
    
    def Msw(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.Msw(t, theta)).T
    
    def Mej(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.Mej(t, theta)).T

    def beta_gamma(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.beta_gamma(t, theta)).T
    
    def beta_th(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.beta_th(t, theta)).T
    
    def R(self, t, theta):
        t, theta = np.meshgrid(t, theta)
        return np.array(self.interpolator.R(t, theta)).T

    # ---------- Radiation Related ---------- #
    # flux density [mJy]
    def FluxDensity(self, t, nu, para, model="sync", xtol=0.0, rtol=1e-2):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.FluxDensity(t, nu, xtol, rtol)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False) / mJy

    # linear polarization
    def Pi_lin(self, t, nu, para, model="sync", xtol=0.0, rtol=1e-2):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.Pi_lin(t, nu, xtol, rtol)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False)
    
    # apparent superluminal motion [mas]
    def Offset(self, t, nu, para, model="sync", xtol=0.0, rtol=1e-2):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.Offset(t, nu, xtol, rtol)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False)
    
    # x size [mas]
    def Xscale(self, t, nu, para, model="sync", xtol=0.0, rtol=1e-2):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.Xscale(t, nu, xtol, rtol)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False)
    
    # y size [mas]
    def Yscale(self, t, nu, para, model="sync", xtol=0.0, rtol=1e-2):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.Yscale(t, nu, xtol, rtol)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False)

    # intensity [cgs]
    def Intensity(self, t, nu, theta, phi, para, model="sync"):
        para_rad = para
        
        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.intensity(t, nu, theta, phi)
        except Exception as e:
            raise e
        
        return np.array(result, copy=False)
    
    def test_function(self, t, nu, theta, phi, para, model="sync"):
        para_rad = para
        
        # configure radiation
        self.afterglow.config_integrator(model, para_rad)

        try:
            result = self.afterglow.test_function(t, nu, theta, phi)
        except Exception as e:
            raise e
        
        return np.asarray(result)

    # imaging
    def solve_image(self, t, nu, para, inum=100, pnum=30, model="sync"):
        para_rad = para

        # configure radiation
        self.afterglow.config_integrator(model, para_rad)
        self.afterglow.config_image()

        # get image data
        try:
            self.afterglow.solve_image(t, nu, inum, pnum)
        except Exception as e:
            raise e
        image = self.afterglow.image()

        return image
