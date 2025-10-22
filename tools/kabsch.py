import numpy as np
from typing import Optional


def kabsch(r1: np.ndarray, r2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotate r2 about the origin the origin to minimize the RMSD from r1.
    
    :param r1: reference structure as an (N, 3) array, where N is the number of atoms
    :type r1: np.ndarray
    :param r2: structure to be aligned as an (N, 3) array
    :type r2: np.ndarray
    :return: A tuple containing:
             - rotation matrix R (3x3 numpy array) that aligns r2 to r1,
             aligned_r2 = r2 @ R
             - aligned coordinates: r2 rotated by the rotation matrix, as an (N, 3) array
    :rtype: tuple[np.ndarray, np.ndarray]
    """
    
    H = np.dot(r2.T, r1)
    U, S, Vt = np.linalg.svd(A)
    if np.linalg.det(np.dot(U, Vt)) < 0.0:
        U[:, -1] *= -1
        R = np.dot(U, Vt)
    else:
        R = np.dot(U, Vt)
    new_r2 = np.dot(r2, R)
    
    return R, new_r2
    