Transition tube sampling
########################

Contributed by Krystof Brezina on Oct 18, 2024.

How to use TTS
**************

Transition tube sampling (TTS) is a method of sampling thermal reactive candidate geometries for, primarily, the training of reactive neural network potentials (NNPs). [1]_ [2]_
In general, this is quite a difficult task, because in order to overcome a thermally insurmountable potential barrier, one has to perform a kind of enhanced sampling simulation to drive the system over it:
this comes with an inpractically high computational cost when done at an *ab initio* level.
TTS significantly reduces this computational cost as it only requires a discretized estimate of the reaction path (*e.g.*, a minimum energy path, MEP) 
and a small numbers of Hessian matrices evaluated at chosen control points (geometries :math:`\mathbf{R}_c, \ c = 1, \dots, N_c`) along the path.
Typically, selecting control points at the reagent and the product minima and the transition state is enough, but one can go as dense as they want.

Setting up
==========

Having calculated the Hesians for the :math:`N_c` control points (and arranged them in an Iterable, *e.g.*, a list `hessians`), one can easily set up the normal mode objects using the `NormalModes` class of the TTS module::

    import transition_tube_sampling as tts
    normalmodes = [tts.NormalModes(h, masses) for h in hessians]

where `masses` is an Iterable of atomic masses used for mass-weighing purposes.
The module has some clever ways of dealing with Hessians from CP2K and ASE, for any other Hessian one must convert it externally to mass-weighted atomic units and 
enter it directly following the documentation below.

Once the local modes `normalmodes` are set up, one can move to the TTS itself through the `NormalModeSampling` class, which works in the following way. 
Given the geometries of the MEP listed in an iterable (*e.g.*, `positions`), the code runs a cubic spline through it to have an analytic description of the MEP
and its tangent vectors (derivatives). 
The other output to set up the `NormalModeSampling` class if the control hessians from above, just note that the code requires a dictionary, 
where the keys are *integer indices* of the control points within the discretized MEP (something like `idx_c = [0, 5, 10]` for a control point at the 0th, 5th and 10th position of the MEP).
You can get this easily through, for instance, a dictionary comprehension::

    assert len(normalmodes) == len(idx_c)
    dnormalmodes = {idx: nm for idx, nm in zip(idx_c, normalmodes)}

and now, the `NormalModeSampling` object is set up just using::

    sampling = tts.NormalModeSampling(positions, dnormalmodes)

and is ready to sample new geometries.

Sampling new geometries
=======================

This is achieved through calling the main interface of the class::

   new_geometries = sampling.get_samples(...)

with appropriate parameters.
The way this works is that for each of the :math:`c` control points, this erects a distribution of *reference geometries* :math:`{\mathbf{R}_0}` peaking at the position
of the control point and decaying away.
This encourages using local modes "locally"; summing up over these distributions gives an exactly uniform distribution along the path.
Now, each of the reference geometries is displaced perpendicularly from the path using its assigned control normal modes

.. math::
    \mathbf{R} = \mathbf{R}_0 + \mathcal{M}^{-\frac{1}{2}}\sum_{i=1}^{N_\mathrm{vib}} \Omega_{i,c} \mathbb{\Omega}_{i,c},

where :math:`\mathcal{M}` is the diagonal mass matrix and the values of normal coordinates :math:`\Omega_{i,c}` are sampled from the harmonic Boltzmann distribution

.. math::
    p(\Omega_{1,c}, \dots, \Omega_{N_\mathrm{vib},c})
    =
    \prod_{i=1}^{N_\mathrm{vib}} \mathrm{exp} \left( -\frac{1}{2}\beta\omega_{i,c}^2 \Omega_{i,c}^2 \right)

at an inverse temperature :math:`\beta = (k_\mathbf{B} T)^{-1}`.
Note that imaginary-frequency modes are ommited as they correspond to the reactive direction which is already accounted for.
Choosing `sampling_mode='quantum'` invokes the calculation of the effective quantum inverse temperature

.. math::
    \beta^*
    =
    \frac{2}{\hbar\omega} \tanh \left( \frac{\beta\hbar\omega}{2} \right)

for each mode separately, yielding the quantum harmonic distribution instead.

Important notes
===============

TTS reduces to a normal mode sampling around a single minimum if `len(dnormalmodes) == 1`.
This can be used to get thermal structures of non-reactive molecules and crystals as long as they are reasonably described by the harmonic approximation.\
In that case, the obtained set of structures is canonical on the harmonic potential energy surface.

In the case of TTS beyond a single minimum, the obtained ensemble of structures does not correspond to the canonical ensemble.
No expectation values should be calculated over it, its purpose is only auxiliary to provide training structures for NNP generation.

Code documentation
******************

.. automodule:: tts.transition_tube_sampling
   :members: NormalModes, NormalModeSampling

References
**********

.. [1] Schran C., Brezina K., Marsalek O. *Committee neural network potentials control generalization errors and enable active learning*. J. Chem. Phys. 153(10), 104105, **2020**
.. [2] Brezina K., Beck H. Marsalek O. *Reducing the cost of neural network potential generation for reactive molecular systems* J. Chem. Theory Comput. 19(19), 6589-6604, **2023**
