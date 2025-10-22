####################################
Time Correlation Functions
####################################

*author*: Elia Stocco

Many physical observables depend on the (Fourier transform of the) Time Correlation Function (TCF) of some microscopic quantities.
Here are some examples:

    - Infrared, Raman and Sum Frequency Generation spectra 
    - frequency dependent dielectric constant and dielectric susceptibility 
    - Vibrational Density Of States

.. (from the electric-dipole TCF)
.. (still from the electric-dipole TCF)
.. (nuclear velocities TCF)

The way to compute these quantities depends on their definition of course, 
but a general way to evaluate the TCF from the time series of microscopic quantities wil lbe provided in this following.

.. attention::
    This page is not complete.
    Here is a list of things that will be added in the future:
        - windowing
        - padding
        - practical examples and code

************************************
Cross Correlation
************************************
The Time Correlation Functions (TCF) between two observables :math:`A,B` is defined through the `cross correlation <https://en.wikipedia.org/wiki/Cross-correlation>`_ (similar, but not the same of the `convolutuon <https://en.wikipedia.org/wiki/Convolution>`_) of their time series:

.. math::

    (A \star B)(\tau) \triangleq & \int_{-\infty}^{\infty} \overline{A(t)} B(t+\tau) \, dt \\
                                                 \triangleq & \int_{-\infty}^{\infty} \overline{A(t-\tau)} B(t)\,dt
                                         

A naive implementation of the previous formula would have a computation cost scaling as :math:`\mathcal{O}\left(N^2\right)` where :math:`N` is the size of the arrays (assumed to be of the same length for simplicity).

However, we can exploit an analogous of the `convolution theorem <https://en.wikipedia.org/wiki/Convolution_theorem>`_ which allows to express the TCF as:

.. math::
    A \star B = \mathcal{F}^{-1} \left[ \overline{\mathcal{F}}\left[ A \right] \cdot \mathcal{F}\left[ B \right] \right]

where :math:`\mathcal{F}\left[ \cdot \right]` is the Fourier transform of a function defined as:

.. math::
    \mathcal{F}[ f ](\omega) \triangleq \int_{-\infty}^{+\infty} e^{-i \omega t} f(t) \, dt 


| The Fourier transfrom has a computational cost of :math:`\mathcal{O}\left(N\log N\right)` thanks to the `Fast Foruier Transform (FFT) <https://en.wikipedia.org/wiki/Fast_Fourier_transform>`_ implementations.
| For this reason, we will provide a way to evaluate the TCF using this "Fourier transform trick".
| Here is a simple function that returns the TCF of an array with itself:

Function Documentation
-----------------------

.. automodule:: tcf.tcf.autocorrelate
   :noindex:

.. literalinclude:: ../tcf/tcf.py
   :language: python

Usage Examples
-----------------------

.. code-block:: python

    import numpy as np
    from tcf import autocorrelate

    # Create an array with 1000 random numbers
    arr = np.random.rand(1000)

    # Compute the autocorrelation function
    autocorr = autocorrelate(arr)

    # Plot the autocorrelation
    import matplotlib.pyplot as plt
    plt.plot(autocorr)
    plt.xlabel('Lag')
    plt.ylabel('Autocorrelation')
    plt.show()

.. attention::
    In practical calculations, we do not have infinite time series with the time :math:`t` spanning in :math:`(-\infty,+\infty)`, 
    but the time ranges from :math:`t \in [0,T]` only.

    Moreover, the numerical routines actually implement a monolateral Fourier transform :math:`\tilde{\mathcal{F}}[ \cdot ]`
    of the following form:

    .. math::
        \tilde{\mathcal{F}}[ f ](\omega) \triangleq \int_{0}^{T} e^{-i \omega t} f(t) \, dt 


************************************
Infrared Spectrum
************************************
The vibrational infrared absorption spectrum :math:`\mathcal{I}\left(\omega\right)`, at thermal equilibrium, can be evaluated using the time series of the dipole using the following expression:

.. math::
    \mathcal{I}\left(\omega\right) 
        = \alpha\left(\omega\right) n\left(\omega\right) 
        & = \, \frac{\pi \omega}{3 \hbar c \Omega \varepsilon_0} 
        \left( 1 - e^{-\beta \hbar \omega} \right) 
        \, \mathcal{C}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right) \\
        & = \, \frac{\pi \omega^2}{3 \hbar c \Omega \varepsilon_0} 
        \, \tilde{\mathcal{C}}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right)

where 
    - :math:`\alpha\left(\omega\right)` is the Beer-Lambert absorption coefficient, 
    - :math:`n\left(\omega\right)` is the refractive index of the material, 
    - :math:`\beta = \frac{1}{k_B T}` is the thermodynamic beta, 
    - :math:`c` is the speed of light, 
    - :math:`\Omega` is the system volume
    - :math:`\varepsilon_0` is the vacuum permittivity
    - :math:`\mathcal{C}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right)` is the Fourier transform of the standard correlation function :math:`\mathcal{C}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(t\right)` [1]_
    - :math:`\tilde{\mathcal{C}}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right)` is the Fourier transform of the Kubo-transformed correlation function :math:`\tilde{\mathcal{C}}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(t\right)`

and :math:`\tilde{\mathcal{C}}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right)` is approximated by using molecular dynamics with the cross correlation of the dipole with itself:

.. math::
    \tilde{\mathcal{C}}_{\boldsymbol{\mu}\boldsymbol{\mu}}\left(\omega\right) \approx 
    \mathcal{F}\left[\boldsymbol{\mu}\star\boldsymbol{\mu}\right]\left(\omega\right) 
    \triangleq \int_{-\infty}^{+\infty} e^{-i \omega t} \boldsymbol{\mu}\star\boldsymbol{\mu}(t) \, dt 

It is important to stress a couple of things:
    - the TCF should actually be calculated for the fluctuation of the dipole w.r.t. to its average value:

    .. math::
        \boldsymbol{\mu}\star\boldsymbol{\mu}(t) \longrightarrow C_{\delta\boldsymbol{\mu},\delta\boldsymbol{\mu}}(t) 
        \quad \text{with} \quad
        \delta\boldsymbol{\mu} = \boldsymbol{\mu} - \braket{\boldsymbol{\mu}}
    and this means that, before computing any Fourier transform, one should remove the mean value of the time series.
    However, we will keep using the same notation as before for simplicity.

    - the definition of the infrared spectrum requires an integration over time which actually range in :math:`(-\infty,+\infty)`.
    The TCF of the dipole with itself is also time-reversible, i.e.

    .. math::
        \left(\boldsymbol{\mu}\star\boldsymbol{\mu}\right)(t) = \left(\boldsymbol{\mu}\star\boldsymbol{\mu}\right)(-t)
    and this implies that its Fourier transform in purely real.

    However, in numerical implementation we can only evaluate the monolateral Fourier transform :math:`\tilde{\mathcal{F}}[ \cdot ]`, whose output is complex since no integration over time for :math:`t\lt 0` occurs.

    For this reason, to guarantee the correct calculation of the infrared spectrum, what is actually implemented numerically is

    .. math::
        \mathcal{I}\left(\omega\right) = &
        \frac{\beta\omega^2}{6c\Omega\varepsilon_0} 
        \, 2 \, {\rm Re} \tilde{\mathcal{F}}\left[\boldsymbol{\mu}\star\boldsymbol{\mu}\right]\left(\omega\right) \\
        = & 
        \frac{\beta\omega^2}{3c\Omega\varepsilon_0} 
        {\rm Re} 
        \int_{0}^{+\infty} e^{-i \omega t} \left(\boldsymbol{\mu}\star\boldsymbol{\mu}\right)\left(t\right) \, dt 

************************************
A common trick
************************************
| Numerical simulations providing the dynamics over time of the dipole (as well as any other quantity) can be used to evaluate the previous formulas only if the simulat time is long enough and the dynamics is properly captured, i.e. if the integration time step is small enough.
| In certain systems, like liquid water, properly sampling the long time behavior of :math:`\left(\boldsymbol{\mu}\star\boldsymbol{\mu}\right)\left(t\right)`, i.e. the low frequency behavior of its Fourier transform, can be challenging due to its intrisically "low dynamics".
| This means that to properly sample the low frequency region of :math:`\mathcal{I}\left(\omega\right)`, extremely long simulations would be needed.
| There is however a possible, or partial, solution that can be applied to the infrared spectrum, which consists in relating :math:`\left(\dot{\boldsymbol{\mu}}\star\dot{\boldsymbol{\mu}}\right)\left(t\right)` with :math:`\boldsymbol{\mu}\star\boldsymbol{\mu}(t)`, i.e. the TCF of the dipole time derivative.
| It is not hard to show that:

.. math::
    \omega^2 \mathcal{F}\left[\boldsymbol{\mu}\star\boldsymbol{\mu}\right]\left(\omega\right) = 
    \mathcal{F}\left[\dot{\boldsymbol{\mu}}\star\dot{\boldsymbol{\mu}}\right]\left(\omega\right)

.. attention::
    Just for completeness, it is worth mentioning that this relation is valid only in the continuous case, while in the discrete case some care should be taken when evaluating the derivative of a function [2]_. 
    You can find more information online under the topic of "modified wave number".
    However, we are neglecting this aspect here for simplicity.

| This relation shows that the TCF of the dipole time derivative decays faster, and we can also use this TCF to evaluate the infrared spectrum:
    
    .. math::
        \mathcal{I}\left(\omega\right)
        = 
        \frac{\beta}{3c\Omega\varepsilon_0} 
        {\rm Re} 
        \int_{0}^{+\infty} e^{-i \omega t} \left(\dot{\boldsymbol{\mu}}\star\dot{\boldsymbol{\mu}}\right)\left(t\right) \, dt 



.. [1] B. J. Braams, T. F. Miller, and D. E. Manolopoulos, Sum rule constraints on Kubo-transformed correlation functions, Chemical Physics Letters 418, 179 (2006)
.. [2] Sunaina, Butola, M., & Khare, K. (2018). Calculating numerical derivatives using Fourier transform: Some pitfalls and how to avoid them. European Journal of Physics, 39(6), 065806