Customized Radiation Models
===========================
It's also possible to define your own radiation model in `jetsimpy`. To do so, you need to define the solid angle integrated specific intensity (:math:`I_{\nu}`) of each fluid element in the comoving frame. You may find the defination in Rybicki & Lightman. This value can be calculated by :math:`I_\nu=j_\nu \Delta R`, where :math:`j_\nu` is the isotropic emissivity, and :math:`\Delta R` is the blast shell width in the coming frame (defined as ``blast.dR`` in the code).

There are two ways to define the radiation model. The pure Python way and the C++ way. The first way is very convinient. You don't have to learn C++, and you don't need to compile the code everytime after you change your model. However the code will be 10x slower than the second way. In the second way you need to define the intensity by a C++ function, and you will have to install the code every time you make any modifications.

The best practice is to do the Python way to quickly prototype and test your model. Then you can do the C++ way in actual MCMC fitting.

The Python way
--------------
First, you need to define a function describing the specific intensity of a fluid element with the following functional form::

    def custom_intensity(nu, P, blast):
        ...
        return intensity

This function must strictly take three parameters: ``nu`` the frequency, ``P`` the parameter dictionary, and ``blast`` the object to access the fluid properties. 

To apply this model in calculating the flux density, you just need to ::

    jet.FluxDensity(t, nu, P, model=custom_intensity)

An example of customized radiation model can be found in the example folder of the source code.

When defining your radiation model, the specific intensity may depend on the properties of the shock. For example, the cooling frequency of the spectrum depends on the evolution time scale :math:`t`, and the optical depth depends on the shell thickness :math:`\Delta R`. These values can be accessed from the function argument ``blast``. Currently, the following shock properties are supported::

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
The C++ way is fundimentally similar to the Python way, where you need to define a "lambda function" about the specific intensity of the fluid element. You can define your model in the file "jetsimpy/src/Afterglow/models.cpp". In this file there is a function::

    void Models::registerIntensity() {
        ...
    }

where you can **define your radiation model as a lambda function in the function body**. All radiation models are saved in a dictionary ``radiation_models`` with a corresponding keyword. To add your own model, you must define a function with the following form::

    radiation_models["model_keyword"] = [&](const double nu, const Dict& P, const Blast& blast) {
        ...
        return intensity;
    };

where ``"model_keyword"`` is a string specified by the user to identify the model. The arguments of this function are defined in the same way as in the Python way. 

After you have defined your model, you must install the code again (described in the starting page) to make it work. To apply this model in calculating the flux density, you just need to ::

    jet.FluxDensity(t, nu, P, model="model_keyword")

In fact, all built-in models are defined in the same way, where you may find them helpful.