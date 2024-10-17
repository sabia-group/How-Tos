####################################
Normal and Phonon Modes
####################################

A simple explanation of:
    - what normal and phonon modes are
    - how to compute them
    - how to "project" a Molecular Dynamics trajectory onto these modes.

We will first give a mathematical derivation of the normal modes of a molecule.
Later on, we will extend this concept to phonon modes of a crystal with a minimal effort.

************************************
Normal Modes
************************************

======================
The Harmonic Approximation
======================

| Let's consider the system moving around a minimum of the Potential Energy Surface (PES).
| The nuclear coordinates of the system will be indicated with :math:`\mathbf{R}_{eq} = \left( R^1_x, R^1_y, R^1_z, R^2_x, \dots , R^N_z \right)`, which is a vector containing all the degrees of freedom (nuclear coordinates) of the sysyem.

For the sake of simplicity we will indicate the displacement on the nuclei w.r.t. :math:`\mathbf{R}_{eq}` just with :math:`\mathbf{R}` instead of :math:`\Delta \mathbf{R} = \mathbf{R} - \mathbf{R}_{eq}`

Let's express the potential energy of the system within the harmonic approximation around the equilibrium positions :math:`\mathbf{R}_{eq}`.

.. math::

    U \approx U_{\text{harm}} = \mathbf{R}^T \cdot \Phi \cdot \mathbf{R}

.. U \approx U_{\text{harm}} = 
.. \frac{1}{2} \sum_{IJ}^{N_{n}} \sum_{\alpha,\beta=x,y,z} 
.. \Delta \mathbf{R}^I_\alpha \cdot \Phi_{IJ}^{\alpha\beta} \cdot \Delta \mathbf{R}^J_\beta 



where :math:`\Phi` is called "force constants matrix" and is defined as:

.. math::
    \Phi_{ij} = \left. \frac{\partial^2 U}{\partial \mathbf{R}_i \partial \mathbf{R}_j }\right|_{\mathbf{R}=\mathbf{R}_{eq}}

.. \Phi_{IJ}^{\alpha\beta} = \left. \frac{\partial^2 U}{\partial \mathbf{R}^I_\alpha \partial \mathbf{R}^J_\beta}\right|_{\mathbf{R}=\mathbf{R}_{eq}}

.. note::

    The indices :math:`i,j` are multi-indices making the previous equations really compact.
    
    For reference, usually the harmonic potential energy and the force constant matrix are expressed in the following way (equivalent but more complicated way):

    .. math::

        \begin{aligned}
        U_{\text{harm}} & = 
        \frac{1}{2} \sum_{IJ}^{N_{n}} \sum_{\alpha,\beta=x,y,z} 
        \Delta \mathbf{R}^I_\alpha \cdot \Phi_{IJ}^{\alpha\beta} \cdot \Delta \mathbf{R}^J_\beta \\
        \Phi_{IJ}^{\alpha\beta} & = \left. \frac{\partial^2 U}{\partial \mathbf{R}^I_\alpha \partial \mathbf{R}^J_\beta}\right|_{\mathbf{R}=\mathbf{R}_{eq}}
        \end{aligned}

    Where indices :math:`I,J\in[1,\dots,N]` label the nuclei and the indices :math:`\alpha,\beta\in\left\{x,y,z\right\}` label the Cartesian direction.



In the following, we will focus on isolated systems only (Γ-only treatment).

We want to find the solutions of the following linearized/harmonic equations of motion:

.. math::

    \mathcal{H} = \sum_I^{N_n} \frac{m^I \mathbf{v}^{I} \cdot \mathbf{v}^{I}  }{2} + 
    \frac{1}{2} \sum_{IJ}^{N_{n}} \sum_{\alpha,\beta=x,y,z} 
    \Delta \mathbf{R}^I_\alpha \cdot \Phi_{IJ}^{\alpha\beta} \cdot \Delta \mathbf{R}^J_\beta 

.. math::

    m^I


