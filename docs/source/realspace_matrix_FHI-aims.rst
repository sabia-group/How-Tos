####################################
Sparse Real-Space Matrix in FHI-aims
####################################

Overview
--------

FHI-aims stores sparse matrices (Hamiltonian and overlap) using a compressed format similar to CSR (Compressed Sparse Row) due to the excellent locality of the real-space representation. The sparse storage significantly reduces memory usage by only storing non-zero elements.

Basic CSR Format Background
----------------------------

Standard CSR Components
~~~~~~~~~~~~~~~~~~~~~~~

* **V**: Array containing all non-zero elements
* **col_index**: Column indices of non-zero elements
* **row_index**: Array indicating where each row starts in the V array

CSR Example
~~~~~~~~~~~

For a 4×4 matrix with 4 non-zero elements (array index starts from 0 in this section):

.. code-block:: text

   [5 0 0 0]
   [0 8 0 0]
   [0 0 3 0]
   [0 6 0 0]

CSR representation:

* ``V = [5, 8, 3, 6]`` (non-zero values)
* ``col_index = [0, 1, 2, 1]`` (column indices of non-zero elements)
* ``row_index = [0, 1, 2, 3, 4]`` (cumulative count: row i has ``row_index[i+1] - row_index[i]`` elements)

FHI-aims Specific Implementation
---------------------------------

Key Differences from Standard CSR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **PBC Consideration**: Due to periodic boundary conditions, the ``row_index`` definition is slightly modified
2. **3D Index Structure**: Uses 3D arrays to handle cell and basis indexing simultaneously
3. **Hamiltonian-Specific Storage**: ``index_hamiltonian`` array structure

Data Structure Components
~~~~~~~~~~~~~~~~~~~~~~~~~~

For Hamiltonian Matrix Storage
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are the FORTRAN variables used in FHI-aims for storing the Hamiltonian:

* **hamiltonian**: Array containing non-zero matrix elements (equivalent to V in CSR)
* **column_index_hamiltonian**: Column indices of stored elements
* **index_hamiltonian**: 3D array with dimensions ``(2, n_cells, n_basis)``
  * First dimension (2): stores start and end indices for each row
  * Second dimension: cell index
  * Third dimension: basis function index

Index Mapping
~~~~~~~~~~~~~

Cell and Basis Organization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* :math:`m` cells, and :math:`p` basis functions per cell
* Full Hamiltonian would be :math:`(m \times p) \times (m \times p)` dimensions
* Due to PBC, only need to store :math:`m \times p^2` elements (cells are equivalent under translation)

Index Structure Example
^^^^^^^^^^^^^^^^^^^^^^^

Array index starts from 1 in this section to follow FORTRAN standard.

For period=3 (so 3 cells in total), and 2 basis functions per cell (think of graphene but in 1D):

.. code-block:: text

   | b1 b2 | b1 b2 | b1 b2 |  (b -> basis)
   |   c3  |   c1  |   c2  |  (c -> cell)

The variable ``cell_index`` in machine storage is:

.. code-block:: text

   cell_index =
    0  0  0  # cell 1: central
    1  0  0  # cell 2: on the right
   -1  0  0  # cell 3: on the left

The full Hamiltonian is (index format ``(cell, basis)``):

.. code-block:: text

          (1,1) (1,2)
   (1,1)     v1    v2  # index_hamiltonian(:,1,1) = [1,2]
   (1,2)    v2*     0  # index_hamiltonian(:,1,2) = [0,-1]
   (2,1)      0    v3  # index_hamiltonian(:,2,1) = [3,3]
   (2,2)      0    v4  # index_hamiltonian(:,2,2) = [4,4]
   (3,1)      0     0  # index_hamiltonian(:,3,1) = [0,-1]
   (3,2)    v3*   v4'  # index_hamiltonian(:,3,2) = [5,5]

So the variable ``hamiltonian`` in machine storage is:

* ``v2*`` is the complex conjugation of ``v2``
* ``v2*`` and ``v3*`` are ignored because of symmetry, and for ``v4'`` although it is equal to ``v4*``, it should still be stored because of being in the upper triangle part

.. code-block:: text

   hamiltonian = [v1, v2, v3, v4, v4']

And the variable ``column_index_hamiltonian`` in machine storage is (it represents the column index of ``hamiltonian``):

.. code-block:: text

   column_index_hamiltonian = [1, 2, 2, 2, 2]

And the variable ``index_hamiltonian`` is:

.. code-block:: text

   index_hamiltonian(1,:,:) =
           basis_1 basis_2
   cell_1        1       0
   cell_2        3       4
   cell_3        0       5

   index_hamiltonian(2,:,:) =
           basis_1 basis_2
   cell_1        2      -1
   cell_2        3       4
   cell_3       -1       5

The following form makes it easier:

.. code-block:: text

   index_hamiltonian(:,1,1) = [1,2]    # row 1 cell_1 basis_1, elements from index 1 to 2
   index_hamiltonian(:,1,2) = [0,-1]   # row 2 cell_1 basis_2, no elements (empty)
   index_hamiltonian(:,2,1) = [3,3]    # row 3 cell_2 basis_1, element at index 3
   index_hamiltonian(:,2,2) = [4,4]    # row 4 cell_2 basis_2, element at index 4
   index_hamiltonian(:,3,1) = [0,-1]   # row 5 cell_3 basis_1, no elements (empty)
   index_hamiltonian(:,3,2) = [5,5]    # row 6 cell_3 basis_2, element at index 5

Special Conventions
~~~~~~~~~~~~~~~~~~~

Empty Rows
^^^^^^^^^^

* If a row has no matrix elements to store: ``index_hamiltonian(:,m,n) = (0,-1)``
* This indicates an empty row with no non-zero elements

Symmetry Considerations
^^^^^^^^^^^^^^^^^^^^^^^

* In the main cell (0, 0, 0)
  
  * due to symmetry: :math:`\langle c_0 b_i | A | c_0 b_j \rangle = \langle c_0 b_j | A | c_0 b_i \rangle^*`
  * we only store the upper half of the matrix, that's why ``v2*`` is ignored

* For two "conjugated" cell I and J which ``cell_index(I,:) = - cell_index(J,:)``
  
  * because of translation symmetry: :math:`\langle c_0 b_i | A | c_I b_j \rangle = \langle c_J b_i | A | c_0 b_j \rangle`
  * only the upper part, like :math:`\langle c_0 b_i | A | c_I b_j \rangle` where :math:`i<j`
    
    * (so this is U in each sub-matrix between cell 0 and cell :math:`I`)
  
  * and the lower part is the complex conjugation, like the case for ``v3*``
  * The reason for storing ``v4'`` is that we have :math:`i=j` and the diagonal part will be stored since it locates in the upper triangle

Data Access Pattern
^^^^^^^^^^^^^^^^^^^

* To extract elements for a specific row: use the index range specified in ``index_hamiltonian``
* Extract ``hamiltonian[start:end+1]`` and corresponding ``column_index_hamiltonian[start:end+1]``

Special row
^^^^^^^^^^^

* The last row in ``cell_index`` is not used
  
  * actually it is ``999999999  999999999  999999999``

* Just skip it, as you can see in the following example code

Use Case in FORTRAN
~~~~~~~~~~~~~~~~~~~

.. code-block:: fortran

   do i_cell_row = 1, n_cells_in_hamiltonian-1   ! yes, "-1".
       do i_basis_row = 1, n_basis
       i_index_first = index_hamiltonian(1, i_cell_row, i_basis_row)
       i_index_last = index_hamiltonian(2, i_cell_row, i_basis_row)
       do i_index = i_index_first, i_index_last
           i_basis_col = column_index_hamiltonian(i_index)
           ! Use:
           !    hamiltonian(i_index, i_spin)
           !    density_matrix_sparse(i_index)
           !    and i_basis_row, i_cell_row, i_basis_col
           ! or any combination of
           !    (i_basis_row, i_loc_cell_row), &
           !    & (i_basis_col, i_loc_cell_col),
           ! with
           !    i_cell_row == &
           !    & position_in_hamiltonian(i_loc_cell_row, i_loc_cell_col)
       end do
       end do
   end do

Key Points for Parser Implementation
-------------------------------------

1. **Index Convention**: Check if FHI-aims uses 0-based or 1-based indexing
2. **Empty Row Handling**: Properly handle ``(0,-1)`` markers for empty rows
3. **PBC Structure**: Account for the cell-basis double indexing system
4. **Symmetry**: Determine if full or half-storage is used for symmetric matrices

   * AND: upper matrix in FORTRAN is lower matrix in python, because of array indexing sequence

5. **Data Types**: Verify floating-point precision and integer types used
6. **FORTRAN / Python indexing**: FORTRAN starting from 1 and include both start and end, but python start from 0 and only include the start

