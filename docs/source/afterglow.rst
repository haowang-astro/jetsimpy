Afterglow modeling
==================
The GRB afterglow observables can be calculated by calling the methods of the jet object `jetsimpy.Jet`.

Light curve and spectrum
------------------------
.. py:method:: .FluxDensity(t, nu, P, model='sync', rtol=1e-3, max_iter=100, force_return=True)

    Calculate the flux density.

    :param np.array/float t: time series (s)
    :param np.array/float nu: frequency series (Hz)
    :param dict P: parameter dictionary
    :param str model: radiation model
    :param float rtol: relative tolerance
    :param int max_iter: adaptive integration maximum iteration number.
    :param int force_return: return the result without throwing an error even if the adaptive integration doesn't achieve the desired accuracy after the maximum iteration number.
    :return: flux density (mJy)

If `t` is a 1D numpy.array and `nu` as a scalar, a light curve is generated. If `t` is a scalar and `nu` as a 1D numpy.array, a spectrum is generated. `t` and `nu` can also be 1D numpy.array of the same length, which is useful in data fitting where different data points have different frequency.

The parameter dictionary `P` must be compatible with the radiation model (to be explained in the next section). For the default model `model="sync"` the required keyword parameters are `eps_e` (:math:`\epsilon_{\rm e}`), `eps_b` (:math:`\epsilon_{\rm B}`), `p`, `theta_v` (:math:`\theta_{\rm obs}`), `d` (luminosity distance), and `z` (redshift).

Frequency integrated flux
-------------------------
.. py:method:: .Flux(t, nu1, nu2, P, model='sync', rtol=1e-3, max_iter=100, force_return=True)

    Calculate the frequency integrated flux.

    :param np.array/float nu1: frequency integration lower limit (Hz)
    :param np.array/float nu2: frequency integration upper limit (Hz)
    :return: flux (erg/s/cm^2)

Apparent superluminal motion
----------------------------
The flux centroid offset can calculated by the following method

.. py:method:: .Offset(t, nu, P, model='sync', rtol=1e-3, max_iter=100, force_return=True)

    Calculate the flux centroid offset. 

    :return: flux centroid offset (MAS)

Image size
----------
.. py:method:: .SizeX(t, nu, P, model='sync', rtol=1e-3, max_iter=100, force_return=True)

    Calculate the Gaussian equivalent 1-sigma image size along the jet axis. 

    :return: x direction image size (MAS)

.. py:method:: .SizeY(t, nu, P, model='sync', rtol=1e-3, max_iter=100, force_return=True)

    Calculate the Gaussian equivalent 1-sigma image size perpendicular to the jet axis. 

    :return: y direction image size (MAS)

Sky map
-------
.. py:method:: .IntensityOfPixel(t, nu, x_offset, y_offset, P, model="sync")
    
    The intensity of a "pixel" in the image.

    :param float t: time series (s)
    :param float nu: frequency series (Hz)
    :param np.array x_offset: offset of the pixel along jet axis direction (MAS).
    :param np.array y_offset: offset of the pixel perpendicular to the jet axis direction (MAS).
    :param dict P: parameter dictionary
    :param str model: radiation model