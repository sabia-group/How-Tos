################################
Nudged Elastic Band Calculations
################################

Contributed by George Trenins

This page describes how to set up an NEB calculation using MPI to parallelize over images. The energy and force calculations are delegated to "client" codes, such as
``FHI-aims``, while a "driver" code performs the actual optimization given the externally calculated energies and forces.

.. toctree::
    :maxdepth: 1

    neb-with-ase