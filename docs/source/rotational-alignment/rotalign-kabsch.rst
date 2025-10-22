################
Kabsch algorithm
################

Contributed by George Trenins

**************
Method summary
**************

The algorithm can be derived as follow. First, expand the square in the definition of the RMSD

.. math::

    \mathbf{R} = \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\min}} \sum_{\alpha=1}^{N} \left \lvert \mathbf{r}_{1,\alpha} - \mathbf{O} \mathbf{r}_{2,\alpha} \right \rvert^2 = 
     \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\min}}
     \sum_{\alpha=1}^{N} \left\{ \left \lvert \mathbf{r}_{1,\alpha} \right \rvert^2
     + \left \lvert \mathbf{r}_{2,\alpha} \right \rvert^2  - 2 \mathbf{r}_{1,\alpha}^{\top}  \mathbf{O} \mathbf{r}_{2,\alpha} \right\}

Therefore, we can equivalently require


.. math::

    \mathbf{R} = \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\max}}
     \sum_{\alpha=1}^{N} \mathbf{r}_{1,\alpha}^{\top}  \mathbf{O} \mathbf{r}_{2,\alpha}  = 
      \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\max}} \ \mathrm{Tr} \left[  \mathbf{O} \sum_{\alpha=1}^{N} \left\{ \mathbf{r}_{2,\alpha} \mathbf{r}_{1,\alpha}^{\top}  \right\} \right].


Defining the :math:`3\times 3` matrix :math:`\mathbf{H} = \sum_{\alpha=1}^{N}  \mathbf{r}_{2,\alpha} \mathbf{r}_{1,\alpha}^{\top}` we apply SVD, :math:`\mathbf{H} = \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{\top}`. The matrix :math:`\mathbf{\Sigma}` is diagonal, real, and non-negative. From the properties of the trace, we may write

.. math::

    \mathbf{R} =  \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\max}} \ \mathrm{Tr} \left[  \mathbf{O}  \mathbf{U} \mathbf{\Sigma} \mathbf{V}^{\top} \right] = 
     \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\max}} \ \mathrm{Tr} \left[  \mathbf{V}^{\top} \mathbf{O}  \mathbf{U} \mathbf{\Sigma} \right] = 
     \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\max}} \ \mathrm{Tr} \left[  \mathbf{S} \mathbf{\Sigma} \right] 

where :math:`\mathbf{S} = \mathbf{V}^{\top} \mathbf{O}  \mathbf{U}`. This matrix is orthogonal by construction.
It is straightforward to find the optimal rotation in this representation:

  * if :math:`\mathrm{det} (\mathbf{\mathbf{U}} \mathbf{V}^{\top}) = 1` then the optimal rotation satisfies :math:`\mathbf{S} = \mathbf{V}^{\top} \mathbf{R} \mathbf{U} = \mathbf{I}`, so that :math:`\mathbf{R} = \mathbf{V} \mathbf{U}^{\top}`.

  * if :math:`\mathrm{det} (\mathbf{\mathbf{U}} \mathbf{V}^{\top}) = -1` then the optimal rotation satisfies :math:`\mathbf{S} = \mathbf{V}^{\top} \mathbf{R} \mathbf{U} = \mathbf{J}`, with :math:`\mathbf{J} = \mathrm{diag}[1,\, \ 1,\, {-1}]`, so that :math:`\mathbf{R} = \mathbf{V} \mathbf{J} \mathbf{U}^{\top}`.


**********************
Function Documentation
**********************

.. autofunction:: tools.kabsch.kabsch
   :noindex:

.. literalinclude:: ../../../tools/kabsch.py
   :language: python


Usage Examples
--------------

.. code-block:: python

    import numpy as np
    from tools.kabsch import kabsch
    from scipy.spatial.transform import Rotation as R

    rng = np.random.default_rng(31415)
    N = 10   # number of "atoms"
    # Generate random structure
    r1 = rng.uniform(low=-5, high=5, size=(N,3))
    # Random rotation matrix
    rot = R.random(random_state=rng).as_matrix()
    r2 = r1 @ rot
    new_rot, new_r2 = kabsch(r1, r2)
    # These should be the same
    print(f"{r1 = }")
    print(f"{new_r2 = }")
    # ...also these
    print(f"{rot.T = }")
    print(f"{new_rot = }")
