####################################
Time Correlation Functions
####################################

Many physical observables depend on the (Fourier transform of the) Time Correlation Functions (TCF) of some microscopic quantities.
Here are some examples:

    - Infrared, Raman and Sum Frequency Generation spectra 
    - frequency dependent dielectric constant and dielectric susceptibility 
    - Vibrational Density Of States

.. (from the electric-dipole TCF)
.. (still from the electric-dipole TCF)
.. (nuclear velocities TCF)

The way to compute these quantities depends on their definition of course, 
but a general way to evaluate the TCF from the time series of microscopic quantities wil lbe provided in this following.


************************************
Cross correlation
************************************
The Time Correlation Functions (TCF) :math:`C_{AB}\left(t\right)` between two observables :math:`A,B` is defined through the `cross correlation <https://en.wikipedia.org/wiki/Cross-correlation>`_ (similar, but not the same of the `convolutuon <https://en.wikipedia.org/wiki/Convolution>`_) of their time series:

.. math::

    C_{AB}\left(\tau\right) = (A \star B)(\tau) \triangleq & \int_{-\infty}^{\infty} \overline{A(t)} B(t+\tau) \, dt \\
                                                 \triangleq & \int_{-\infty}^{\infty} \overline{A(t-\tau)} B(t)\,dt
                                         

A naive implementation of the previous formula would have a computation cost scaling as :math:`\mathcal{O}\left(N\right)` where :math:`N` is the size of the arrays (assumed to be of the same length for simplicity).

However, we can exploit an analogous of the `convolution theorem <https://en.wikipedia.org/wiki/Convolution_theorem>`_ which allows to express the TCF as:

.. math::
    A \star B = \mathcal{F}^{-1} \left[ \overline{\mathcal{F}}\left[ A \right] \cdot \mathcal{F}\left[ B \right] \right]

where :math:`\mathcal{F}\left[ \cdot \right]` is the Fourier transform of a function defined as:

.. math::
    \mathcal{F}\left[ f \right](\omega) \triangleq \int_{-\infty}^{+\infty} e^{-i \omega t} f\left(t\right) \, dt

The Fourier transfrom has a computational cost of :math:`\mathcal{O}\left(N\log N\right)` thanks to the `Fast Foruier Transform (FFT) <https://en.wikipedia.org/wiki/Fast_Fourier_transform>`_ implementations.

For this reason, we will provide a way to evaluate the TCF using this "Fourier transform trick".