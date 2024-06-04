Customized Radiation Models
===========================
It's also possible to define your own radiation model in `jetsimpy`. To do so, you need to define the emissivity (:math:`\epsilon'_{\nu'}`) of each fluid element in the comoving frame. 

There are two ways to define the emissivity function. The pure Python way and the C++ way. The first way is very convinient. You don't have to learn C++, and you don't need to compile the code everytime after you change your model. However the code will be 10x slower than the second way. In the second way you need to define the emissivity by a C++ function, and the speed is quite fast.

The best practice is to do the Python way to quickly prototype and test your model, and do the C++ way in actual MCMC fitting.

The Python way
--------------
First, you need to define an emissivity function with the following function form::

    def custom_emissivity(nu, P, blast):
        ...
        return emissivity

The function must strictly take three parameters: ``nu`` the frequency, ``P`` the parameter dictionary, and ``blast`` an object to access the fluid properties. 

Then, you just need to set the keyword argument ``model`` to the function name::

    jet.FluxDensity(t, nu, P, model=custom_emissivity)

An example of customized emissivity can be found in the example folder of the source code.

The ``blast`` object has the following properties::

    # coordinate values (burster frame)
    blast.t                 # time since burst
    blast.theta             # polar angle of the fluid element
    blast.phi               # azimuthal angle of the fluid element
    blast.R                 # Radius of the fluid element

    # blast velocity (burster frame)
    blast.beta             # (post-shock) velocity
    blast.gamma            # (post-shock) Lorentz factor
    blast.beta_th          # (post-shock) polar velocity component
    blast.beta_r           # (post-shock) radial velocity component
    blast.beta_f           # (forward shock) velocity
    blast.gamma_f          # (forward shock) Lorentz factor
    blast.s                # calibration coefficient
    blast.doppler          # Doppler factor

    # thermodynamic values (comoving frame)
    blast.n_blast          # (post-shock) number density
    blast.e_density        # (post-shock) energy density (rest mass excluded)
    blast.pressure         # (post-shock) pressure
    blast.dR               # shell width

    # others (burster frame)
    blast.n_ambient        # external number density

The C++ way
-----------
