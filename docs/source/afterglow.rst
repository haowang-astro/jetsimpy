Afterglow modeling
==================
The GRB afterglow observables can be calculated by calling the methods of the jet object `jetsimpy.Jet`.

Light curve and spectrum
------------------------
.. py:method:: .FluxDensity(t, nu, P, model='sync', rtol=1e-3)

    Calculate the flux density.

    :param np.array/float t: time series (s)
    :param np.array/float nu: frequency series (Hz)
    :param dict P: parameter dictionary
    :param str model: emissivity model
    :param float rtol: relative tolerance
    :return: flux density (mJy)

If `t` is a 1D numpy.array and `nu` as a scalar, a light curve is generated. If `t` is a scalar and `nu` as a 1D numpy.array, a spectrum is generated. `t` and `nu` can also be 1D numpy.array of the same length, which is useful in data fitting where different data points have different frequency.

The parameter dictionary `P` must be compatible with the emissivity model (to be explained in x). For the default model `model="sync"` the required keyword parameters are `eps_e` (:math:`\epsilon_{\rm e}`), `eps_b` (:math:`\epsilon_{\rm B}`), `p`, `theta_v` (:math:`\theta_{\rm obs}`), `d` (luminosity distance), and `z` (redshift).

Apparent superluminal motion
----------------------------
The flux centroid offset can calculated by following method

.. py:method:: .Offset(t, nu, P, model='sync', rtol=1e-3)

    Calculate the flux centroid offset. 

    :return: flux centroid offset (MAS)

Image size
----------
.. py:method:: .SizeX(t, nu, P, model='sync', rtol=1e-3)

    Calculate the Gaussian equivalent 1-sigma image size along the jet axis. 

    :return: x direction image size (MAS)

.. py:method:: .SizeY(t, nu, P, model='sync', rtol=1e-3)

    Calculate the Gaussian equivalent 1-sigma image size perpendicular to the jet axis. 

    :return: y direction image size (MAS)

Sky map
-------
.. py:method:: .IntensityOfPixel(t, nu, x_tilde, y_tilde, P, model="sync")
    
    The intensity of a "pixel" which is offset from the burst center by (x_tilde, y_tilde).