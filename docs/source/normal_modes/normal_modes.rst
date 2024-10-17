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

| For the sake of simplicity we will indicate the displacement on the nuclei w.r.t. :math:`\mathbf{R}_{eq}` just with :math:`\mathbf{R}` instead of :math:`\Delta \mathbf{R} = \mathbf{R} - \mathbf{R}_{eq}`

| Let's express the potential energy of the system within the harmonic approximation around the equilibrium positions :math:`\mathbf{R}_{eq}`.

.. math::

    U \approx U_{\rm harm } = \frac{1}{2} \mathbf{R}^T \cdot \Phi \cdot \mathbf{R}

.. U \approx U_{\text{harm}} = 
.. \frac{1}{2} \sum_{IJ}^{N_{n}} \sum_{\alpha,\beta=x,y,z} 
.. \Delta \mathbf{R}^I_\alpha \cdot \Phi_{IJ}^{\alpha\beta} \cdot \Delta \mathbf{R}^J_\beta 

| where :math:`\Phi` is called "force constants matrix" and is defined as:

.. math::
    \Phi_{ij} = \left. \frac{\partial^2 U}{\partial \mathbf{R}_i \partial \mathbf{R}_j }\right|_{\mathbf{R}=\mathbf{R}_{eq}}

.. \Phi_{IJ}^{\alpha\beta} = \left. \frac{\partial^2 U}{\partial \mathbf{R}^I_\alpha \partial \mathbf{R}^J_\beta}\right|_{\mathbf{R}=\mathbf{R}_{eq}}

| where the indices :math:`i,j` runs over all the degrees of freedoms, i.e. :math:`i,j\in[1,\dots,3N]`.
| This means that :math:`\Phi` is a :math:`3N\times3N` matrix.
.. note::

    .. The indices :math:`i,j` can be understood as multi-indices making the previous equations really compact.
    
    | For reference, usually the harmonic potential energy and the force constant matrix are expressed in the following way (equivalent but more complicated way):

    .. math::

        \begin{aligned}
        U_{\text{harm}} & = 
        \frac{1}{2} \sum_{IJ}^{N_{n}} \sum_{\alpha,\beta=x,y,z} 
        \Delta \mathbf{R}^I_\alpha \cdot \Phi_{IJ}^{\alpha\beta} \cdot \Delta \mathbf{R}^J_\beta \\
        \Phi_{IJ}^{\alpha\beta} & = \left. \frac{\partial^2 U}{\partial \mathbf{R}^I_\alpha \partial \mathbf{R}^J_\beta}\right|_{\mathbf{R}=\mathbf{R}_{eq}}
        \end{aligned}

    | where the indices :math:`I,J\in[1,\dots,N]` label the nuclei and the indices :math:`\alpha,\beta\in\left\{x,y,z\right\}` label the Cartesian direction.
    | However, the previous expression using the indices :math:`i,j` are equivalent but much more compact.


We want to find the solutions of the (Newton) equations of motion of following linearized/harmonic hamiltonian:

.. math::

    \mathcal{H} = & T + U_{\rm harm } \\
                = & \frac{1}{2}   \mathbf{v}^T \cdot \mathbf{M} \cdot \mathbf{v} + \frac{1}{2} \mathbf{R}^T \cdot \Phi \cdot \mathbf{R}
    
    
where the kinetic energy :math:`T` has been expressed using the same convention previously adopted for the displacement and the potential energy,
and :math:`\mathbf{M}` is a **diagonal** :math:`3N\times3N` matrix containing the masses of all the atoms (repeated 3 times):

.. math::

    {\rm diag} \, \mathbf{M} = \left( M^1, M^1, M^1, M^2, \dots , M^N \right)


