##########
ASE driver
##########

Contributed by George Trenins

.. warning::
    
    The procedure outlined below was only tested for ``ase-3.25.0``.

The python script for driving the simulation and all accompanying input files can be downloaded :download:`here <../_static/neb_ase_aims.tar.gz>`. Below I explain the key parts of the ``main()`` section in the python script.

**********
Band setup
**********

The nudged elastic band (NEB) method locates the minimum-energy path between two configurations. We assume that these have already been generated and can be read from the files :file:`reactant.in` and :file:`product.in`. The initial setup is standard and follows `ASE documentation <https://wiki.fysik.dtu.dk/ase/ase/neb.html>`_.

.. code:: python

    from ase.io import read
    from ase import Atoms
    from ase.mep import NEB

    reactant = "reactant.in"
    product = "product.in"
    num_images = 4   # number of images between the end-points 
    initial: Atoms = read(reactant, format="aims")
    final: Atoms = read(product, format="aims")
    configs: list[Atoms] = [initial] + [initial.copy() for i in range(num_images)] + [final]
    band: NEB = NEB(configs, parallel=True, method="improvedtangent")
    band.interpolate(method='linear', mic=True, apply_constraint=True)

****************
Socket interface
****************

By far the most efficient way for ``ase`` and ``FHI-aims`` to communicate is via `sockets <https://wiki.fysik.dtu.dk/ase/ase/calculators/socketio/socketio.html>`_. Setting this up requires some care, especially if using multiple nodes.


Driver hostname
===============

In what follows, we will need to advertise the actual hostname of the driver node (node running ``ase``) to the clients. When using the Slurm workload manager, this can be obtained as follows:

.. code:: python
    
    import os, subprocess
    nodelist = os.environ['SLURM_NODELIST']
    nodes = subprocess.check_output(
        ['scontrol', 'show', 'hostnames', nodelist],
        text=True
    ).split()
    driver_host = nodes[0]

.. hint::

    If everything is running on the same node you may simply use :code:`driver_host = "localhost"`.


Launcher scripts
================

There are several ways to configure how the calculator gets launched. I recommend manually instantiating a "profile":

.. code:: python

    from ase.calculators.aims import Aims, AimsProfile
    profile = AimsProfile("/path/to/launcher_script.sh")
    aims = Aims(profile=profile, ...) 

Here, :file:`launcher_script.sh` should be an executable file that launches the calculator for the image (bead) to which it is attached. We will generate the launcher scripts from within the python script driving the calculation and tailor their contents image by image, as I explain below.


Port assignment
===============

For a typical NEB optimization, I recommend using TCP/IP sockets (as opposed to UNIX sockets), since the communication overhead is certain to be negligible compared to the time needed for the *ab initio* calculation, and since this is the more flexible option. To open such a socket, you need to assign an integer `port number <https://en.wikipedia.org/wiki/List_of_TCP_and_UDP_port_numbers>`_ that is not in use by any other applications. See the :code:`get_port()` function implemented in the FHI-aims software package in  :file:`utilities/get_free_port.py` for an example of how to generate a suitable port number.

Calculator initialization
=========================

The exterior images (reactant and product) are treated differently to the interior images in NEB optimization, so we consider the two separately.

Reactant and product
^^^^^^^^^^^^^^^^^^^^

These structures do not change over the course of the optimization. For the most basic NEB optimization method (:code:`method = "aseneb"`) these images do not need a calculator attached at all. All other methods require the potential energies, but not the forces. Since the structures do not change, it is sufficient to compute the energies once and then close the calculator, best accomplished using a context manager. The computed energies are cached and persist after the calculator is closed. At this stage, 
the contents of :file:`launcher_script.sh`  for these images can be something like

.. code::

    srun /path/to/aims.VERSION.scalapack.mpi.x < /dev/null > aims.out

allowing ``aims`` to utilise all the available CPUs, since we compute the energies first for the reactant and then for the product, so competition for resources is not an issue.

.. code:: python

    for i in [0, num_images + 1]:
        image: Atoms = band.images[i]
        target: Path = Path(f"image{i:02d}")
        port = get_port(host = driver_host)  
        cmd = command_for_exterior_images      # e.g., 'srun /path/to/aims.VERSION.scalapack.mpi.x < /dev/null > aims.out'
        launcher = wd / f"_launcher{i:02d}.sh" # separate launcher for every images
        write_launcher(launcher, cmd)          # see below
        profile = AimsProfile(str(launcher))
        aims = Aims(
            profile=profile,
            directory=target,
            ...,                                # species_dir, kgrid, etc.
            use_pimd_wrapper=(driver_host, port),
        )
        # use context manager to open and close socket
        with SocketIOCalculator(aims, log=sys.stdout, port=port) as calc:
            image.calc = calc
            # no need to assign the energy, cached by image
            image.get_potential_energies()


Interior images
^^^^^^^^^^^^^^^

The energies and forces on the interior images are recomputed at every step of the NEB optimization, therefore we need to:

  #. Launch several calculators in parallel and keep them running until the optimization is done.

  #. Ensure that the available resources are evenly distributed.

The last point in particular has presented some unexpected difficulties (at least on ADA). In my tests using 2 nodes (72 cores each), when launching four MPI processes (36 tasks each), slurm would routinely assign three processes to one node, and only one process to the other. The following did not fix the issue:

  * waiting a few seconds between launching different client processes

  * using the :code:`slurm --exact --exclusive` flag combination 

  * using :code:`slurm --distribution=cyclic`

A robust approach is to manually assign nodes to the different images. The example I give below assumes that the calculator for a single image runs on only one node (not necessarily utilising all the CPUs), but can be readily extended to multi-node jobs. In this case, the launcher command is something like

.. code::

    srun -N 1 -n 36 --exact --exclusive /path/to/aims.VERSION.scalapack.mpi.x < /dev/null > aims.out

and the python script goes as follows:

.. code:: python

    from itertools import cycle
    import time
    node_cycle = cycle(nodes)  # `nodes` defined in the 'Driver hostname' section
    for i,(image,node) in enumerate(zip(band.images[:-1], node_cycle)):
        if i == 0: continue    # reactant already taken care of
        target: Path = Path(f"image{i:02d}")
        launcher = wd / f"_launcher{i:02d}.sh"
        port = get_port(host = driver_host)
        cmd = command_for_exterior_images  # e.g., srun -N 1 -n 36 --exact --exclusive /path/to/aims.VERSION.scalapack.mpi.x < /dev/null > aims.out
        write_launcher(launcher, cmd, extra=f"--nodelist={node}")  # force round-robin
        profile = AimsProfile(str(launcher))
        aims = Aims(...)  # same as before
        calc = SocketIOCalculator(calc=aims, port=port, log=sys.stdout)
        # Manually launch the client and the server so that get_port() 
        # has up-to-date info on what ports are available in the next 
        # iteration of the for loop
        calc.server = calc.launch_server()
        proc = calc.launch_client(image, properties=["energy", "forces"],
                                  port=calc._port,
                                  unixsocket=calc._unixsocket)
        time.sleep(1.0) # optional
        calc.server.proc = proc 
        image.calc = calc

The function :code:`write_launcher()` used in this and preceding section is

.. code:: python

    from typing import Optional
    import stat 
    def write_launcher(
            filepath: Path, 
            cmd: str, 
            extra: Optional[str] = None) -> None:
        if extra is not None:
            cmd_lst: list[str] = cmd.split()
            cmd_lst.insert(1, extra)
            cmd: str = ' '.join(cmd_lst)
        with open(filepath, 'w') as f:
            f.write(f"#!/bin/bash\n\n{cmd}\n")
        # Get current permissions
        mode = filepath.stat().st_mode
        # Add execute permission for user, group, and others
        filepath.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        return
