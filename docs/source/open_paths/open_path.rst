###########################################
Momentum distributions from open path PIMD
###########################################

Contributed by Hannah Bertschi

This tutorial shows how to set up an open path PIMD simulation in i-pi for the example of the water dimer with an open path on one of the hydrogens. Furthermore, it describes how the output from i-pi has to be used and processed to get the radial momentum distribution. Potential problems and how to test the results are provided as well. 

**********
Background
**********

This section provides very minimal information on what equations are necessary to calculate the momentum distribution of an atom in a multi-atom molecule in three dimensions. For derivations please refer to [Pro]_ the supporting information also contain important information.

The momentum distribution of atom :math:`k` can be represented as a Fourier transform of the off-diagonal position matrix element of the density matrix :math:`e^{- \beta \hat{H}}/Z` 

.. math::
    n_k(\boldsymbol{p}) = \frac{1}{(2 \pi \hbar)^3 Z} \int \ d\boldsymbol{\Delta} \  e^{\frac{i}{\hbar}\boldsymbol{p}\cdot \boldsymbol{\Delta}} \int d\{\boldsymbol{q}\} \bra{\boldsymbol{q}_1, \ldots \boldsymbol{q}_N} \ e^{- \beta \hat{H}} \ket{\boldsymbol{q}_1, \ldots \boldsymbol{q}_k + \boldsymbol{\Delta}, \ldots  \boldsymbol{q}_N}

The subscripts refer to the :math:`N` atoms and :math:`\boldsymbol{q}` is their three dimensional position vector. The matrix element is off-diagonal, as there is a shift for atom :math:`k` of :math:`\boldsymbol{\Delta}`. The position matrix element can be sampled by running path-integral molecular dynamics, where on atom :math:`k` no spring between the first and last bead is present. 

In the end we want a radial momentum distribution, since during the simulation the molecule either way rotates. After averaging over the angles we get 

.. math::
    n_k(p) = 4 \pi \int d\Delta \ \mathcal{N}_k(\Delta) \ \Delta^2 \frac{\hbar}{p \Delta} \sin \left( \frac{p \Delta}{\hbar}\right)

where :math:`\mathcal{N}_k(\Delta)` can be sampled by running an open path PIMD simulation. The distance of the first and last bead :math:`P` of atom :math:`k` given by :math:`|\boldsymbol{q}_k^{(P)} - \boldsymbol{q}_k^{(1)}|` has to be calculated for each sample and then averaged via the following function

.. math::
    \mathcal{N}_k(\Delta) = \biggl \langle \frac{(2 \pi \sigma_P^2)^{-1/2}}{\Delta \ |\boldsymbol{q}_k^{(P)} - \boldsymbol{q}_k^{(1)}|} \left[  e^{-\frac{1}{2 \sigma_P^2} (\Delta - |\boldsymbol{q}_k^{(P)} - \boldsymbol{q}_k^{(1)}|)^2} - e^{-\frac{1}{2 \sigma_P^2} (\Delta + |\boldsymbol{q}_k^{(P)} - \boldsymbol{q}_k^{(1)}|)^2}\right] \biggr \rangle.

Here we use :math:`\sigma_P^2 = \hbar^2 \beta/P m_k` with the inverse temperature :math:`\beta`, the amount of beads :math:`P` and the mass of the atom :math:`m_k`.

************
i-pi inputs
************

Similar to standard PIMD we run a molecular dynamics simulation in the NVT ensemble with multiple beads. Most of the input keywords stay the same as compared to the standard PIMD. We will go here through an example xml input file and highlight what input is special for the open path simulations. 

Some thing to think about before/ when setting up the input file are (*italic* for what is shown in this example):

- What is the temperature of interest? *60 K* 
- How many beads are necessary for the chosen temperature and system? *64 beads*
- What potential to use? *connection to MBX forcefield via unix socket*
- Is there already a PIMD simulation from which to restart? *use RESTART-64 for initialization*
- For which atom do I want to calculate the momentum distribution, i.e. which atom(s) need an open path?
  *one open path on atom 1 (hydrogen)*
- What thermostat and related parameters to use?
- ...


.. note::
   The counting of indices in i-pi starts from zero. This is relevant for specifying the atom with an open path.

Lets first look at what output to request:

.. code-block:: xml

  <simulation verbosity='low'>
    <output prefix="simulation">
      <properties stride='40' filename='out'>
            [ step, time{picosecond}, conserved, temperature{kelvin}, kinetic_cv, potential, kinetic_cv(2), kinetic_cv(1) ]
      </properties>
      <properties stride='40' filename='H.q'> [ atom_x_path(1)] </properties>
      <trajectory bead='-16' format='xyz' filename='pos' stride='40'> positions </trajectory>
      <trajectory bead='63' format='xyz' filename='posn' stride='40'> positions </trajectory>
      <checkpoint filename="chk" stride="2000" overwrite="true"/>
    </output>

The momentum distribution depends on the distance between the first and last bead of the atom with the open path. Therefore these geometries need to be ouptput. This is done in two ways here, the line ``<properties stride='40' filename='H.q'> [ atom_x_path(1)] </properties>``. It outputs all the bead positions of atom 1 every 40 steps.  

Other than writing all bead positions of the atom with open path we can also write all atoms of some specific beads. The line ``<trajectory bead='-16' format='xyz' filename='pos' stride='40'> positions </trajectory>`` outputs every 16nth bead of each atom, i.e. including bead 0. And the next line outputs the geometry of the last bead 63. Each bead is output in a xyz file.

Either option gives us the geometries of atom 1 for the first and last bead. These we will need to process in the following section to get the momentum distributions.

Next in the xml file follows some generic information on the total steps and connection to a driver.

.. code-block:: xml

  <total_steps>100000</total_steps>
    <prng>
      <seed>3348</seed>
    </prng>
    <ffsocket mode='unix' name='driver'>
      <address>mbx</address>
    </ffsocket>

Here we tell i-pi to use 64 beads and read for initialization a restart file, which corresponds to an equilibrated ring polymer structure. The forces are just the ones from the driver.

.. code-block:: xml

  <system>
    <initialize nbeads='64'>
      <file mode='chk'> RESTART-64 </file>
    </initialize>
    <forces>
      <force forcefield='driver'/>
    </forces>

This code block specifies the open path. In this case it is on atom 1 (counting from zero), which is a hydrgen atom.

.. code-block:: xml

  <normal_modes>
    <open_paths> [1] </open_paths>
  </normal_modes>

Lastly follows all information on the NVT ensemble. 

.. code-block:: xml

    <ensemble>
      <temperature units='kelvin'>60.0</temperature>
    </ensemble>
    <motion mode='dynamics'>
      <dynamics mode='nvt' splitting='baoab'>
        <thermostat mode='pile_l'>
          <tau units='femtosecond'> 100 </tau>
        </thermostat>
        <timestep units='femtosecond'>0.25</timestep>
      </dynamics>
    </motion>
    </system>
    </simulation>

********************************
Processing the simulation output
********************************

**References**

.. [Pro] V. Kapil, A. Cuzzocrea, and M. Ceriotti. *Anisotropy of the Proton Momentum Distribution in Water* J. Phys. Chem. B **122** 6048-6054 (2018).

