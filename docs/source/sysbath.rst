Dissipative Baths
=================

Basics
------

A dissipative environment can be represented by a set of harmonic oscillators, which for a one-dimensional system with phase-space coordinates :math:`(P,\,Q)` takes the from

.. math::

    H(P,Q,p,q) = \frac{P^2}{2 m} + V(Q) + \sum_{n = 1}^{N_{\text{bath}}} 
    \frac{p_n^2}{2 m} +
    \frac{c_n^2}{2 m \omega_n^2} \left(
        q_n - \frac{c_n Q}{m \omega_n^2}
    \right)^2

The bath frequencies :math:`\omega_n` and coupling coefficients :math:`c_n` come from the discretization of the spectral density :math:`J(\omega)`

.. math::

    J(\omega) \sim \frac{\pi}{2} \sum_{n=1}^{N_{\text{bath}}} 
    \frac{c_n^2}{m \omega_n} \delta(\omega - \omega_n)

which we have implemented for several common models.

Bath discretization
-------------------

.. automodule:: tools.baths
    
    .. autoclass:: tools.baths.BaseSpectralDensity
        :class-doc-from: init
        :member-order: bysource
        :members:

    .. autoclass:: tools.baths.ExpOhmic
        :class-doc-from: init
        :show-inheritance:


