Hydrodynamic Simulation
=======================
The hydrodynamic evolution of a jet with a given angular profile can be solved by creating an object `jetsimpy.Jet`. Below is 

simulation
----------
.. py:class:: jetsimpy.Jet(profile, nwind, nism, tmin=10, tmax=3.2e9, grid=jetsimpy.Uniform(257), tail=True, spread=True, cal_level=1, rtol=1e-6, cfl=0.9)
    
    Create an object containing the simulation results.

    :param tuple profile: a tuple of the tabulated data of the jet angular profile taking the form of `(theta, Eiso, lf)`. `theta` is the polar angle array starting from 0 and end with pi. `Eiso` is the isotropic equivalent energy. `lf` is the Lorentz factor. All three elements must be 1D `np.array` of the same length.
    :param float nwind: wind density scale.
    :param float nism: ism density scale.
    :param float tmin: start time of the simulation.
    :param float tmax: end time of the simulation.
    :param np.array grid: the polar angles of the "cell edges" of the simulation. This parameter specifies the resolution of the simulation.
    :param bool tail: Add an low-energy isotropic tail or not. This tail ensures that the energy is positive and the Lorentz factor is above unity. 
    :param bool spread: turning on the spreading or not.
    :param int cal_level: calibration level. "`cal_level=0`": no calibration. "`cal_level=1`": calibrate with Blandford-McKee solutions throughout the simulation. "`cal_level=2`": calibrate with Blandford-McKee solution in the relativistic limit and Sedov-Taylor solution in the Newtonian limit. Setting `cal_level=2` might be aggressive since the thin shell approximation is no longer valid in the Newtonian limit.
    :param float rtol: relative tolerance in the primitive variable solver. **Never change this parameter unless you know what is going on.**
    :param float cfl: the CFL number. **Never change this parameter unless you know what is going on.**

jet profile
-----------
The jet profile can be arbitrary but must be physically reasonable. For general purpose, users should always provide tabulated data even if the profile is analytic. To save time, the code provides some built-in profiles.

.. py:function:: jetsimpy.TopHat(theta_c, Eiso, lf0=1e100)

    A top-hat jet profile.

    :param float theta_c: half-opening angle (rad).
    :param float Eiso: isotropic equivalent energy.
    :param lf0: initial Lorentz factor.
    :return: a tuple (theta, Eiso, lf)

.. py:function:: jetsimpy.Gaussian(theta_c, Eiso, lf0=1e100)

    A Gaussian jet profile defined by

    .. math::
        E(\theta) &= E_{\rm iso} \exp\left[{-\frac{1}{2}\left(\frac{\theta}{\theta_{\rm c}}\right)^2}\right]

        \Gamma(\theta) &= (\Gamma_0 - 1)\exp\left[{-\frac{1}{2}\left(\frac{\theta}{\theta_{\rm c}}\right)^2}\right] + 1
    
    :param float theta_c: half-opening angle (:math:`\theta_{\rm c}`).
    :param float Eiso: isotropic equivalent energy (:math:`E_{\rm iso}`).
    :param float lf0: initial Lorentz factor (:math:`\Gamma_0`).
    :return: a tuple (theta, Eiso, lf)

.. py:function:: jetsimpy.PowerLaw(theta_c, Eiso, lf0=1e100, s=4.0)

    A power-law jet profile defined by

    .. math:: 
        E(\theta) &= E_{\rm iso} \left[1 + \left(\frac{\theta}{\theta_{\rm c}}\right)^2\right]^{-s/2}

        \Gamma(\theta) &=(\Gamma_0 - 1)\left[1 + \left(\frac{\theta}{\theta_{\rm c}}\right)^2\right]^{-s/2} + 1

    :param float theta_c: half-opening angle (:math:`\theta_{\rm c}`).
    :param float Eiso: isotropic equivalent energy (:math:`E_{\rm iso}`).
    :param float lf0: initial Lorentz factor (:math:`\Gamma_0`).
    :param float s: slope (:math:`s`).
    :return: a tuple (theta, Eiso, lf)

coasting phase
--------------
To turn off the coasting phase, just set a very large (but finite) initial Lorentz factor such as :math:`10^{100}`.

external density
----------------
The external density is assumed to follow the profile below.

.. math::
    n = n_{\rm ism} + n_{\rm wind}(r/10^{17}{\rm cm})^{-2}

start and end time
------------------
The start time `tmin` is non-zero because the initial radius must be positive. The code will expand the blastwave at the initial speed until `tmin`, and then start the simulation. In principle `tmin` can be arbitrarily small. However, to speed up the simulation, a moderate value is preferred. `tmin=10` is in general a safe value.

The simulation ends at `tmax`. This value should be sufficiently large to cover the time span of afterglow modeling.

resolution
----------
In a finite volume method, the simulation domain is discretized into a number of cells. In this code the discretization is determined by the `grid` parameter. The `grid` parameter is a 1D numpy array of the "cell edges". Namely, the cell centers locate at the middle of adjacent elements. 

The resolution directly determines the speed of the code. The optimal way to choose the resolution is to set a high resolution inside the jet core and low resolution outside the core. However, to maintain numerical stability, the size of adjacent cells shouldn't change dramatically. 

The code provides some built-in resolution setup for jet profiles with well-defined opening angles.

.. py:function:: jetsimpy.Uniform(theta_c, npoints)

    Uniform cell distribution. A conservative choice for all cases.

    :param float theta_c: half-opening angle
    :param int npoints: number of cell edges
    :return: np.array: array of cell edges

.. py:function:: jetsimpy.NorthPole(theta_c, npoints)

    Resolution optimized for forward-jet only. 

    :param float theta_c: half-opening angle
    :param int npoints: number of cell edges
    :return: np.array: array of cell edges

.. py:function:: jetsimpy.SouthPole(theta_c, npoints)

    Resolution optimized for counter-jet only. 

    :param float theta_c: half-opening angle
    :param int npoints: number of cell edges
    :return: np.array: array of cell edges

.. py:function:: jetsimpy.BothPoles(theta_c, npoints)

    Resolution optimized for forward-jet and counter-jet.

    :param float theta_c: half-opening angle
    :param int npoints: number of cell edges
    :return: np.array: array of cell edges

simulation results
------------------
One can access the simulation data from the `jetsimpy.Jet` object.

.. py:property:: .t_pde
    
    time data (1d numpy.array)

.. py:property:: .y_pde

    primitive variable data with the shape (5, ntheta, ntime). The 5 variables are [:math:`M_{\rm sw}`, :math:`M_{\rm ej}`, :math:`\beta^2\gamma^2`, :math:`\beta_{\theta}`, :math:`R`].

.. py:property:: .theta_pde

    array of cell centers

.. py:method:: .dE0_dOmega(t, theta)

    The interpolated energy profile (rest mass excluded).

    :param numpy.array t: time
    :param numpy.array theta: polar angle

.. py:method:: .dMsw_dOmega(t, theta)

    The interpolated swept-up mass profile.

    :param numpy.array t: time
    :param numpy.array theta: polar angle

.. py:method:: .dMej_dOmega(t, theta)

    The interpolated ejecta mass profile.

    :param numpy.array t: time
    :param numpy.array theta: polar angle

.. py:method:: .beta_gamma(t, theta)

    The interpolated 4-velocity profile.

    :param numpy.array t: time
    :param numpy.array theta: polar angle

.. py:method:: .beta_theta(t, theta)

    The interpolated tangent velocity profile.

    :param numpy.array t: time
    :param numpy.array theta: polar angle

.. py:method:: .R(t, theta)

    The interpolated Radius angular profile.

    :param numpy.array t: time
    :param numpy.array theta: polar angle