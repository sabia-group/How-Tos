Fluctuations
=====

During a molecular dynamics simulation one generally focuses on three main ensembles, namely 
the :math:`NVT`, :math:`NPT` and :math:`NVE` ensembles.
 
In the case of :math:`NVT` one introduces a thermostat to **equilibrate** the temperaturexi, whereas in :math:`NPT` a barostat is also
used. Various methods for introducing both a thermostat as well as a barostat have been introduced.  

It is important to point out, however, that the target temperature and/or pressure 

.. _Ideal Gas:

Let us consider an ideal gas consisting of :math:`N` particles, then 

:math:`PV = Nk_B T`.

From a classical point of view, if the Hamiltonian is quadratic with :math:`\nu` degress of freedom

:math:`E = \frac{\nu}{2}k_B T`

We emmidiately get that

:math:`E = \frac{2}{3} PV`

.. _note :: 

The proportionality relation is, in fact, more general, and can be applied to poliatomic systems as well.

.. _installation:

Installation
------------

To use Lumache, first install it using pip:

.. code-block:: console

   (.venv) $ pip install lumache

Creating recipes
----------------

To retrieve a list of random ingredients,
you can use the ``lumache.get_random_ingredients()`` function:

.. autofunction:: lumache.get_random_ingredients

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']

