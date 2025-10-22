###########################################
Momentum distributions from open path PIMD
###########################################

Contributed by Hannah Bertschi

This tutorial shows how to set up an open path PIMD simulation in i-pi for the example of the water dimer with an open path on one of the hydrogens. Furthermore, it describes how the output from i-pi has to be used and processed to get the radial momentum distribution. Potential problems and how to test the results are provided as well. 

************
i-pi inputs
************

Similar to standard PIMD we run a molecular dynamics simulation in the NVT ensemble with multiple beads. Most of the input keywords stay the same as compared to the standard PIMD. We will go here through an example xml input file and highlight what input is special for the open path simulations. 

Some thing to think about before/ when setting up the input file are (italic for what is shown in this example):

- What is the temperature of interest? *60 K* 
- How many beads are necessary for the chosen temperature and system? *64 beads*
- What potential to use? *connection to MBX forcefield via unix socket*
- Is there already a PIMD simulation from which to restart? *use RESTART-64 for initialization*
- For which atom do I want to calculate the momentum distribution, i.e. which atom(s) need an open path?
  *one open path on atom 1 (hydrogen)*
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

Other than writing all bead positions of the atom with open path we can also write all atoms of some specific beads. The line ``<trajectory bead='-16' format='xyz' filename='pos' stride='40'> positions </trajectory>`` outputs every 16nth bead of each atom, i.e. including bead 0. And the next line outputs the geometry of the last bead 63.

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



