#################################
Rotational Alignment of Molecules
#################################

Contributed by George Trenins

This tutorial explains how to determine the optimal rotation matrix :math:`\mathbf{R}` that minimizes the root-mean square distance between two molecular structures,


.. math::

    \mathbf{R} = \mathrm{\arg}\,\underset{\mathbf{O}}{\mathrm{\min}} \sum_{\alpha=1}^{N} \left \lvert \mathbf{r}_{1,\alpha} - \mathbf{O} \mathbf{r}_{2,\alpha} \right \rvert^2


There are (at least) two ways of doing this efficiently:

.. toctree::
    :maxdepth: 1

    rotalign-kabsch
    rotalign-quaternion


.. hint::

    In some cases, e.g., for imposing the rotational Eckart condition, you may be interested in minimizing the mass-weighted RMSD :math:`\sum_{\alpha=1}^{N} m_{\alpha} \left \lvert \mathbf{r}_{1,\alpha} - \mathbf{O} \mathbf{r}_{2,\alpha} \right \rvert^2`. The same algorithms can be used for this if the input positions :math:`\mathbf{r}_{i,\alpha}` are replaced with mass-weighted positions :math:`\mathbf{q}_{i,\alpha} = \sqrt{m_{\alpha}} \mathbf{r}_{i,\alpha}`.

.. warning::
    The linked codes do not themselves perform translational alignment. It is the responsibility of the user to 
    align the structures' centroids or centres of mass before feeding them to the rotational alignment code.


