###############################
Creating AFM and STM images
###############################

This is a short workflow regarding the creation of STM and AFM images with different
quantities outputted in :code:`.cube` files from the FHI-aims code. Scripts and
older documents that aid visualization are also linked here.

************************
STM simulations
************************

FHI-aims can easily produce the necessary ingredients to simulate
STM images within the `Tersoff-Hamann approximation`_. Within the framework of
DFT and under some applied voltage bias :math:`V`, the tunneling current shows the following proportionality

.. math::

   I_t(\boldsymbol{r}) \propto
   \begin{cases}
       \sum_{\boldsymbol{k}, \nu} |\psi_{\boldsymbol{k},\nu}(\boldsymbol{r})|^2 f(\varepsilon_{\boldsymbol{k}, \nu} - E_f) [1-f(\varepsilon_{\boldsymbol{k}, \nu} - (E_f-eV))]  & \text{if } V < 0 \\
       \sum_{\boldsymbol{k}, \nu} |\psi_{\boldsymbol{k},\nu}(\boldsymbol{r})|^2 [1-f(\varepsilon_{\boldsymbol{k}, \nu} - E_f)] f(\varepsilon_{\boldsymbol{k}, \nu} - (E_f-eV))  & \text{if } V > 0
   \end{cases}

where :math:`\psi_{\boldsymbol{k}, \nu}` are the Kohn-Sham orbitals with respective energies :math:`\varepsilon_{\boldsymbol{k}, \nu}` and :math:`E_f` is the Fermi energy.
The condition on the occupation functions make sure that the sums only include states within an energy window between :math:`E_f` and :math:`V`.
Note that if :math:`V>0`, the current probes unoccupied states. :math:`\sum_{\boldsymbol{k}, \nu}^{\text{chosen states}} |\psi_{\boldsymbol{k},\nu}(\boldsymbol{r})|^2` is usually called the local
density of states LDOS.

As explained in the manual of FHI-aims, once one has a system geometry for
which an STM image should be created, a full SCF convergence including
the :code:`control.in` command needs to be run:

.. code-block::

  output cube stm 0.3

Note that the last value is given in eV and it can be negative or positive.
It corresponds to the energy window, starting from the Fermi energy determined in the
calculation. At the end of this calculation, a cube file will be generated
where the LDOS is given at each space point in the grid.

There are two possible modes of producing STM images: Constant current and constant height.

For the constant height image, different values of current will be present at each position, and building
an STM image corresponds to plotting all values of the cube file that lie at a given
distance :math:`d` above the sample of interest.
The color code is the value of the LDOS at each point with a fixed `d`.

For the constant current image, building
an STM image corresponds to finding an isosurface in the cube file that maintains the value of the
LDOS constant, within a certain threshold.
The color code is the height of each point for a fixed value of the LDOS.

Because building an image in the constant current case is a bit more involved, but also more common
to be of interest, we provide a script and a few steps to make such images below.

************************
AFM simulations
************************


.. _Tersoff-Hamann approximation: https://link.aps.org/doi/10.1103/PhysRevB.31.805
