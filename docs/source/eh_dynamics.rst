.. _eh-dynamics:
##############################
Ehrenfest dynamics in FHI-aims
##############################

Contributed by Hannah Bertschi

The aim here is to show how to set up an input file in order to run real-time time-dependent density functional theory (RT-TDDFT) calculations together with Ehrenfest dynamics in FHI-aims. This is shown here when an external electric field is applied. The output of the dipole moment could be used to calculate an absorption spectrum, as is shown in this tutorial :ref:`abs-spectrum`. As an example we will look here at the benzene moleucle.

*******************************
Set up the ``control.in`` file
*******************************
As with any density functional theory calculation the functional has to be specified and charges, spin, relativistic corrections, scf parameters etc. can be specified as well in the ``control.in`` file. Here we will use a pbe functional with no spin, charge and relativistic corrections ::

        xc		pz-lda
        spin		none
        relativistic	none
        charge		0.

The propagation of the electronic degrees of freedom is done with RT-TDDFT. Working in atomic units the trajectory is here evaluated for 5500 a.u. with a time step of 0.05 a.u. and output is generated every 0.5 a.u.. The propagation scheme used is ``crank_nicolson`` (first order). ::
        
        RT_TDDFT_input_units atomic
        RT_TDDFT_propagation 5500 0.05 0.5
        RT_TDDFT_propagator crank_nicolson

The propagation of the nuclei with Ehrenfest dynamics is done in the default setting (only one available right now) with the same time steps as the electronic one [#f1]_. ::

        RT_TDDFT_ehrenfest default 0.05 0.5

A time-dependent delta kick pulse is defined like this in FHI-aims

.. math::
   \boldsymbol{E}(t) = \boldsymbol{E}_0 \exp \left(\frac{t - t_0}{t_w}\right) \left( 1 + \exp\left(\frac{t - t_0}{t_w}\right) \right)^{-2}

The field is defined like this ::

        RT_TDDFT_td_field t_start t_end type freq cycle center width Ex Ey Ez

``freq`` and ``cycle`` are not relevant for the delta kick. It has the ``type`` 2, ``center`` corresponds to :math:`t_0` and ``width`` to :math:`t_w`. The three ``Ex``, ``Ey`` and ``Ez`` values make up the :math:`\boldsymbol{E}_0` vector. 

Here we choose ::

        RT_TDDFT_td_field 0 0 2 0 0 10 0.8 0.04 0 0
        RT_TDDFT_td_field_gauge length

This means the pulse is evaluated from the beginning of the trajectory to the end (the first two zeros) and the delta kick pulse is centered around 10 a.u., has a width of 0.8 a.u. and is in x-direction with amplitude 0.04 a.u.. It is applied in the length gauge, which is recommended for molecules.

We can also control what is output to standard output (first letter) and to a separate file (second letter) ::

        RT_TDDFT_output_energies T T
        RT_TDDFT_output_dipole T T
        RT_TDDFT_ehrenfest_output_trajectory T T
        RT_TDDFT_output_level 2


********************************
Set up the ``geometry.in`` file
********************************

.. rubic:: Footnotes
.. [#f1] However outputting only every 0.5 a.u. does not work right now, the Ehrenfest trajectory is output every 0.05 a.u.
