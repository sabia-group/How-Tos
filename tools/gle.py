#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   stochastic.py
@Time    :   2024/06/07 11:46:34
@Author  :   George Trenins
@Contact :   gstrenin@gmail.com
@Desc    :   Some useful functions for stochastic dynamics
'''


from __future__ import print_function, division, absolute_import
import numpy as np
from typing import Union, Optional



def gle_cxx(
        t: Union[float, np.ndarray], 
        beta: float,
        omega: float,
        A_pp: Union[float, np.ndarray],
        mass: Optional[float] = 1.0) -> np.ndarray:
    """NOTE: when calculating the TCF for a WNLE thermostat, simply pass
    gamma = 1/tau for the value of A_pp.

    Args:
        t (Union[float, np.ndarray]): times for which to compute the TCFs
        beta (float): reciprocal temperature
        omega (float): harmonic frequency
        A_pp (Union[float, np.ndarray]): drift matrix
        mass (Optional[float], optional): particle mass. Defaults to 1.0.

    Returns:
        tcf (np.ndarray): 
            shape = (nt, nf, nf) where nt = len(flatten(t)) and nf = len(A_pp)+1
                * tcf[:,0,0] gives cxx
                * tcf[:,1,1] gives cpp
                * tcf[:,0,1] gives cxp
                * etc...
    """
    
    from scipy.linalg import expm
    abs_t = np.reshape(np.abs(t), (-1,1,1))
    A_pp = np.atleast_2d(A_pp)
    if A_pp.ndim != 2:
        raise RuntimeError(f"2-dimensional array expected for A_pp, instead {A_pp.ndim = }")
    if (n := A_pp.shape[0]) != A_pp.shape[1]:
        raise RuntimeError(f"Square A_pp expected, instead got {A_pp.shape = }")
    A_qp = np.zeros((n+1, n+1))
    A_qp[1:,1:] = A_pp
    A_qp[0, 1] = -1
    A_qp[1, 0] = omega**2
    C = np.eye(n+1) / beta
    C[0,0] /= omega**2
    tcf = expm(-A_qp * abs_t) @ C
    sqm = np.sqrt(mass)
    sqmvec = np.ones(n+1)
    sqmvec[0] = 1/sqm
    sqmvec[1] = sqm
    sqmmat = sqmvec[:,None] * sqmvec[None,:]
    return tcf * sqmmat

def wnle_cqq(
        t: Union[float, np.ndarray], 
        beta: float,
        omega: float,
        tau: float,
        m: Optional[float] = 1):
    Omega = np.sqrt(np.abs(omega**2 - 1/(2*tau)**2))
    eps = 1.0e-4
    test = omega*tau - 0.5
    if test > eps:
        ans = np.exp(-t/(2*tau)) * (np.cos(Omega*t) + np.sin(Omega*t) / (2*Omega*tau))
    elif test < -eps:
        from tools.special import logcosh, logsinh
        ans = np.exp(-t/(2*tau) + logcosh(Omega*t))
        ans += np.exp(-t/(2*tau) + logsinh(Omega*t).real) / (2*Omega*tau)
    else:
        ans = np.exp(-t/(2*tau)) * (1 + t/(2*tau))
    return ans / (beta * m * omega**2)