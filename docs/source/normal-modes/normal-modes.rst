Phonons
=======

Author: Elia Stocco

Projection onto Phonon Modes
----------------------------

This section describes the procedure to project molecular dynamics trajectories onto the phonon modes of the system, providing a mode-resolved picture of the dynamics.

The procedure is based on the standard phonon theory of crystals, which can be found in many condensed matter physics textbooks [Rigamonti2007]_ [Ashcroft1978]_. Only the aspects essential for understanding the phonon projection procedure are briefly described here. First, the theory of normal/vibrational modes of a molecule (a system without translational symmetry) is introduced. The extension to crystalline systems and phonons follows.

The Harmonic Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For convenience and practical purposes, nuclear displacement **q** and velocities **v** relative to the equilibrium configuration are written in a compact, "flattened" notation:
  
.. math::
    \mathbf{q} = (q^1_x, q^1_y, q^1_z, q^2_x, \dots, q^{N_a}_x, q^{N_a}_y, q^{N_a}_z)
  
.. math::
    \mathbf{v} = (v^1_x, v^1_y, v^1_z, v^2_x, \dots, v^{N_a}_x, v^{N_a}_y, v^{N_a}_z)

Using the harmonic approximation, the potential energy **U** is expressed as a quadratic function of both displacement **q** and velocities **v**:

.. math::
    \mathcal{H} = \frac{1}{2} \mathbf{v}^\mathrm{t} \cdot \mathbf{M} \cdot \mathbf{v} + \frac{1}{2} \mathbf{q}^\mathrm{t} \cdot \mathbf{\Phi} \cdot \mathbf{q}

where:

.. math::
    \Phi_{ij} = \left. \frac{\partial^2 U}{\partial \mathbf{q}_i \partial \mathbf{q}_j } \right\vert_\text{eq}

Here, **Φ** is the force constants matrix (a square \(3N_a \times 3N_a\) matrix), and **M** is a diagonal \(3N_a \times 3N_a\) matrix containing the atomic masses on the diagonal:

.. math::
    \text{diag}(\mathbf{M}) = (m_1, m_1, m_1, m_2, \dots, m_{N_a}, m_{N_a}, m_{N_a})

The Newton's equations of motion for this Hamiltonian are:

.. math::
    \mathbf{M} \cdot \ddot{\mathbf{q}}(t) = -\mathbf{\Phi} \cdot \mathbf{q}(t)

Assuming an ansatz for coordinates **q**(t):

.. math::
    \mathbf{q}(t) = \mathrm{Re} \, A_n \mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n \cos(\omega_n t)

.. math::
    \mathbf{v}(t) = -\mathrm{Re} \, \omega_n A_n \mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n \sin(\omega_n t)

where **A_n** is a constant with units of \( \mathrm{mg}^{1/2} \). Substituting this ansatz into the equation of motion results in the following eigenvalue problem:

.. math::
    \mathbf{D} \cdot \boldsymbol{\varepsilon} = \boldsymbol{\varepsilon} \cdot \mathbf{\Lambda}

where:

.. math::
    \mathbf{D} = \mathbf{M}^{-1/2} \cdot \mathbf{\Phi} \cdot \mathbf{M}^{-1/2}

is the dynamical matrix, and:

.. math::
    \Lambda_{nm} = \delta_{nm} \omega_n^2

are the eigenfrequencies.

The eigenvectors of **D** are contained in the orthogonal matrix **ε**, and the matrix **Λ** contains the eigenvalues (all non-negative).

The normal modes **N_n** of the system are defined by removing the mass scaling:

.. math::
    \mathbf{N}_n = \frac{\mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n}{\left| \mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n \right|}

It can be shown that the normal modes are normalized but not orthogonal.

Projection onto Normal Modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The exact solution of the molecular dynamics trajectory can still be expressed in terms of the normal modes, even in the anharmonic case:

.. math::
    \mathbf{q}(t) = \sum_n^{3N_a} \mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n \, \tilde{q}_n(t)

.. math::
    \mathbf{v}(t) = \sum_n^{3N_a} \mathbf{M}^{-1/2} \cdot \boldsymbol{\varepsilon}_n \, \omega_n \tilde{v}_n(t)

The coefficients \( \tilde{q}_n(t) \) and \( \tilde{v}_n(t) \) can be calculated by inverting these equations:

.. math::
    \tilde{\mathbf{q}}(t) = \boldsymbol{\varepsilon}^\mathrm{t} \cdot \mathbf{M}^{1/2} \cdot \mathbf{q}(t)

.. math::
    \tilde{\mathbf{v}}(t) = \mathbf{\Lambda}^{-1/2} \cdot \boldsymbol{\varepsilon}^\mathrm{t} \cdot \mathbf{M}^{1/2} \cdot \mathbf{v}(t)

For numerical stability, modes with \( \omega \approx 0 \) (rigid translations/rotations) can be discarded. From these coefficients, one can calculate the energy contributions from each vibrational mode.

References
----------

.. [Rigamonti2007] Rigamonti, A. *Structure and Dynamics* (2007).

.. [Ashcroft1978] Ashcroft, N. W., and Mermin, N. D. *Solid State Physics* (1978).
