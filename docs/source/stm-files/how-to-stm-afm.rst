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
STM images within the Tersoff-Hamann [`approximation`_]. As explained
in the manual of the code, once a one has a system geometry for 
which an STM image should be created, a full SCF convergence including
the :code:`control.in` command needs to be run:

.. code-block::

  output cube stm 0.3

Note that the last value is given in eV, can assume negative or positive values,
 and it corresponds to the energy window, starting from the Fermi energy determined in the 
calculation, over which the electronic local density of states will be integrated. We will
call this quantity the ILDOS. At the end of this calculation, a cube file will be generated
where the ILDOS is given at each space point in the grid. 

There are two possible modes of producing STM images: Constant current and constant voltage.


.. _approximation: https://link.aps.org/doi/10.1103/PhysRevB.31.805
