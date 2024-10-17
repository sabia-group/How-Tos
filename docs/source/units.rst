Unit systems
============


.. automodule:: mdtools.units
    :member-order: bysource
    :members: atomic, hartAng, kcalAfs, kcalAamu, eVAamu

    .. autoclass:: mdtools.units.SI
        :exclude-members:

Physical constants
------------------

You can get the values of some physical constants by instantiating a unit system and accessing the corresponding property

.. code-block:: python

    import units

    atomic = units.atomic()
    print(atomic.hbar)       # 1.0
    print(atomic.amu)        # 1822.8884862173131
    print(atomic.kb)         # 3.1668115634438576e-06

The following constants are implemented in all unit-system classes

.. autoclass:: mdtools.units.SI
    :noindex:
    :exclude-members:
    :member-order: bysource
    :members: hbar, e, kb, amu, me, c


Dimensions
----------

There is a similar mechanism for getting the base dimensional unit converted to SI (implemented dimensions are `energy, length, time, mass, charge, luminous_intensity, amount, current, temperature, action, angular_momentum, vacuum_permitivity`)


.. code-block:: python

    print(atomic.length)     # 5.29177210903e-11
    print(atomic.mass)       # 9.1093837015e-31
    print(atomic.action)     # 1.0545718176461567e-34


Unit conversion
---------------

A "human-readable" unit conversion can be done using the methods

    .. autoclass:: mdtools.units.SI
        :noindex:
        :exclude-members:
        :member-order: bysource
        :members: str2base, str2SI


.. code-block:: python

    atomic.str2base("1 mp")  # mass of proton in atomic units, prints 1836.15267344
    atomic.str2SI("1 a0")    # Bohr radius in metres, prints 5.29177210903e-11
    atomic.str2base("1 fs")  # femtosecond in atomic units of time, 41.341373335335184


Conversion to/from wavenumbers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conversion between base units of energy or radial frequency and wavenumbers (:math:`\text{cm}^{-1}`) is handled separately.

    .. autoclass:: mdtools.units.SI
        :noindex:
        :exclude-members:
        :member-order: bysource
        :members: energy2wn, wn2energy, omega2wn, wn2omega