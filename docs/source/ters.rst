Simulating tip-enhanced Raman spectroscopy (TERS) images
########################################################

Contributed by Krystof Brezina on Oct 22, 2025.

Theory fundamentals
*******************

Tip-enhanced Raman spectroscopy (TERS) is a method of vibrational spectrscopy that relies on the plasmonic reponse of atomically sharp metallic tips to enhance the Raman scattering from single molecules on conductive surfaces.
It can be used for the imaging of vibrational motion by scanning the tip across the surface-bound nanostructures and recording the Raman intensity as a function of the position of the tip: this is known as TERS imaging.
The exact, first-principles  calculation of TERS images is a difficult, time-dependent problem not applicable to systems with realistic extents.
This stems primarily from the need to represent the spatially non-homogeneous, varying plasmonic near fields induced by the incoming laser radiation at the tip-molecule junction.
To enable TERS in realistic system, the group devised an approximate scheme [1]_ that maps the original, complicated problem onto a much simpler one where the plasmonic near field is represented
as a static embedding potential within a standard density-functional theory (DFT) single-point evaluation with the following, perturbed Hamiltonian

.. math::
    \hat{H} = \hat{H}_0 + E_z \left[ -\hat{\mu}_z + \left( \frac{\partial \hat{\Phi}}{\partial E_z} \right)_0 \right].

In this expression, :math:`\hat{H}_0` is the unperturbed Hamiltonian of the surface-molecule system, :math:`E_z` is the electric far-field (laser) component perpendicular to the surface, :math:`\hat{mu}_z` is the corresponding component of the dipole moment operator of the surface-molecule system and, finally, :math:`\hat{\Phi}` is the electrostatic potential of the isolated tip.
All the details and derivations can be found in References [1]_ and [2]_. 
The important point here is that :math:`\hat{\Phi}` only depends on the spatial coordinate and parametrically on the position of the tip relative to the studied system: :math:`\Phi(\mathbf{r}; \mathbf{R})` and, as such, can be used as a simple embedding potential inside a perturbation DFT calculation with a homogeneous external field :math:`E_z`.
Within our infrastructure, it comes numerically tabulated as a Gaussian cube file and is used in an off-the-shelf manner for simulating TERS images of arbitrary surface-molecule systems.
A self-consistent solution of such Hamiltonian yields an electronic density :math:`\rho(\mathbf{r}; \mathbf{R})`, which can be spatially integrated over the unit cell to yield the perpendicular dipole moment component 

.. math::
    \mu_z(\mathbf{R}) = \int_\text{unit cell} \mathrm{d}\mathbf{r} \ z \rho(\mathbf{r}; \mathbf{R}),

which we need for the calculation of the relevant :math:`zz`-component of the Raman tensor.
Additionally, this allows us to use periodic boundary conditions and represent extended metallic surfaces in the :math:`xy`-plane [2]_.
These dipoles are calculated for two values of :math:`E_z` (typically one of them is 0) and a polarizability is obtained using a finite-difference calculation:

.. math::
    \alpha_{zz}(\mathbf{R}) = \frac{\partial \mu_z(\mathbf{R})}{\partial E_z} \approx \frac{\mu_z(E_z, \mathbf{R}) - \mu_z(0, \mathbf{R})}{E_z}.

Finally, the Raman intensity is obtained by performing a finite difference of the polarizability with respect to a normal-mode coordinate :math:`Q_k`:

.. math::
    I_{zz}(\mathbf{R}) = \left( \frac{\partial \alpha_{zz}(\mathbf{R})}{\partial Q_k} \right)^2.

All in all, to calculate a TERS singal for a given tip position :math:`\mathbf{R}` one needs to do 4 single-point evaluations to numerically calculate the mixed second-derivative of :math:`\mu_z(\mathbf{R})` with respect to :math:`E_z` and :math:`Q_k`.


Implementation
**************

The calculation of TERS spectra requires two components.
First, a tool that is capable of performing a single point DFT calculation with a non-homogeneous embedding potential.
We implemented this in FHI-aims.
Second, the FHI-aims calculation must be wrapped into a larger, overarching infrastructure that manages the calculation of the finite differences, normal modes and tip positions.
I implemented this as a Python package available under our group `GitHub <https://github.com/sabia-group/periodic-ters>`__.
Both aspects are desribed below. 

Single-point calculation with embedding near fields
===================================================

A recent (>=06/2025) version of `FHI-aims <https://fhi-aims.org>`__ is needed.
Enabling the near-field embedding as described above is easy and requires new specific keywords in both the `control.in` and `geometry.in` files.
The control file needs the following:

.. code::

    output	 	    dipole
    pos_tip_origin          -0.000030 -1.696604 -4.614000
    pos_sys_origin          0.000000 0.000000 0.000000
    tip_molecule_distance   4.000000
    rel_shift_from_tip      -4.500000 -4.500000
    nearfield_derivative    tipA_05_vh_ft_0049_3221meV_x1000.cube

The `pos_tip_origin` keyword is a fixed parameter (3D position) that marks the position of the tip apex within the provided cube file. This will always be the same for a given Cube file.
The `pos_sys_origin` keyword is a parameter (3D position) that marks the reference origin in the surface-molecule system (typically this would be the center of mass of the molecule). 
A good starting point is to align the system to that position explicitly, FHI-aims uses the Wigner-Seitz cell anyway.
The `tip_molecule_distance` parameter is a 1D length and encodes the distance between the tip apex and `pos_sys_origin`.
The `rel_shift_from_tip` is a lateral 2D displacement in the :math:`xy`-plane that tell us how the tip should be moved relative to the chosen `pos_sys_origin`.
All distance units are in Ångström.
Finally, the `nearfield_derivative` keyword is a string with the name of the Gaussian Cube file that contains the tabulated near field. 
You can download one from our `GitHub repo <https://github.com/sabia-group/periodic-ters>`__ (to be made public soon, tell us if you want access before that happens).
It is necessary to also include the `output dipole` option to have the dipole moment printed out in the output.

The geometry file needs the following

.. code::

    homogeneous_field 0.0 0.0 -<value>

to turn on the external far field.
Note that the minus sign is necessary to due unusual internal sign conventions in FHI-aims electrostatics. This setting triggers an external field of +<value> V/Å as understood in the usual fashion.
Additionally, `lattice_vector` keywords can be added to close the system into periodic boundary conditions. This will only work for systems with surface slabs and will be checked during the calculation.

.. warning::

    If you are dealing with a periodic system, make sure that your unit cell is larger than the Cube file containing the embedding near-field potential.
 

Simulating TERS line spectra (1D) and images (2D)
=================================================

The above fully defines a working FHI-aims calculation, but is impractical once a large number of such calculations is needed.
To make it useful, I prepared the Finite-Field TERS package that takes care of preparing the numerous FHI-aims inputs with all the parameters correctly set and running them in a tidy manner.
This runs all the necessary calculations from a simple Python script and is fully documented below.
Start by adding the code to your Python path (for example by running `source env.sh` on the provided environment file) and importing the object as

.. code::

    import FiniteFieldTERS as ffters

First, one initiates the `FiniteFieldTERS` object with all the information that the calculation needs from the user, such as the Hessian of the system, the value of the external field, the required displacement along the normal modes, *etc*. 
It can look something like this 

.. code::

   ters = ffters(
      hessian = hessian,
      modes = None,
      masses = masses,
      dq = 5e-3,
      efield = -1e-1,
      submit_style = 'slurm',
      fn_control_template = Path('template.in'),
      aims_dir = Path('/u/brek/build/FHIaims/'),
      fn_batch = Path('run.sbatch'),
      fn_tip_derivative = Path('tipA_05_vh_ft_0049_3221meV_x1000.cube'),
      fn_tip_groundstate = None,
      fn_geometry = Path('geometry.in'),
      ) 

The `template.in` is a text-file stub of the FHI-aims input that contains all the usual choices of the xc functional, van-der-Waals corrections, *etc*.
This will be decorated by the keywords introduced above.
Note the `submit_style` keyword: setting it to `draft` only creates all the neat calculation directories, setting it to `slurm` should correctly submit all calculations too under a SLURM batch system at HPC machines, provided a `fn_batch` SLURM run script.
Do let us know is something is broken, please. 

Then, we are typically interested in two regimes of TERS.
In the first one, we want to calculate a 1D (line) spectrum with a fixed tip position over the frequency range. 
This can be accomplished with the `run_1d_multimode()` functionality with a particular call looking like, *e.g.*,

.. code::
    
    ters.run_1d_multimode(
        mode_indices=np.arange(len(ters.modes)),
        tip_origin = (-0.000030, -1.696604, -4.6140),
        sys_origin = (0.0, 0.0, 0.0),
        tip_height = 4.0,
        xy_displacement = (0.0, 0.0),
        dump_wavenumbers=True
        ) 


Note that we are submitting a set of different mode indices and a single `xy_displacement`, which corresponds to a static tip.

In the second regime (and probably the more common one), we aim to calculate a 2D TERS image by moving the tip across the molecule and recording the TERS intensity along the way.
This calculation is requested by calling `run_2d_grid()`, which can look like

.. code::

    ters.run_2d_grid(
        idx_mode = 14,
        tip_origin = (-0.000030, -1.696604, -4.6140),
        sys_origin = (0.0, 0.0, 0.0),
        tip_height = 4.0,
        scan_range = (-5.0, 5.0, -5.0, 5.0),
        bins = (10, 10)
        )

Note that in this case, we provide a single normal mode, but a `scan_range` and `bins` to set up a grid for the various tip positions.

Once the calculations are finished, the package provides analyses functions to gather the data, calculate the required finite differences and assemble the TERS spectra and images in what essentially is a one-liner context.
For inspiration check out the `examples` directory in our `GitHub <https://github.com/sabia-group/periodic-ters>`__: these should give you a practical handle on how to run things.

Code documentation
******************

.. automodule:: ters.finite_field_ters
   :members: FiniteFieldTERS, analyze_1d_ters, analyze_2d_ters

References
**********

.. [1] Litman Y. et al. *First-Principles Simulations of Tip Enhanced Raman Scattering Reveal Active Role of Substrate on High-Resolution Images*. The Journal of Physical Chemistry Letters 14(30), **2023**
.. [2] Brezina K., Litman Y., Rossi M. *Explaining Principles of Tip-Enhanced Raman Images with Ab Initio Modeling*. arXiv:2509.13075, **2025**
