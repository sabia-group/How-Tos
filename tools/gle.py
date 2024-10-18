#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Some analytical results for a harmonic oscillator attahced to a generalised Langevin equation (GLE) thermostat. Written by George Trenins and Hannah Bertschi.
'''


from __future__ import print_function, division, absolute_import
import numpy as np
from typing import Union, Optional


def wnle_cqq(
        t: Union[float, np.ndarray], 
        beta: float,
        omega: float,
        tau: float,
        m: Optional[float] = 1):
    """Position auto-correlation function for a harmonic oscillator at thermal equilibrium, under the action of a white-noise Langevin equation thermostat (WNLE).

    :param t: time(s) for which to compute the covariance matrix
    :type t: numpy.ndarray
    :param beta: reciprocal temperature :math:`1/k_B T`
    :type beta: float, optional
    :param omega: frequency of the harmonic oscillator
    :type omega: float
    :param tau: reciprocal of the friction, :math:`\tau = 1/\gamma`
    :type tau: float
    :param mass: mass of the harmonic oscillator
    :type mass: float, optional
    :return: covariance matrix
    :rtype: numpy.ndarray
    """
    
    Omega = np.sqrt(np.abs(omega**2 - 1/(2*tau)**2))
    eps = 1.0e-4
    test = omega*tau - 0.5
    if test > eps:
        ans = np.exp(-t/(2*tau)) * (np.cos(Omega*t) + np.sin(Omega*t) / (2*Omega*tau))
    elif test < eps:
        from tools.special import logcosh, logsinh
        ans = np.exp(-t/(2*tau) + logcosh(Omega*t))
        ans += np.exp(-t/(2*tau) + logsinh(Omega*t).real) / (2*Omega*tau)
    else:
        ans = np.exp(-t/(2*tau)) * (1 + t/(2*tau))
    return ans / (beta * m * omega**2)


def gle_cxx(
        t: Union[float, np.ndarray], 
        omega: float,
        Ap: Union[float, np.ndarray],
        Cp: Optional[Union[float, np.ndarray]] = None,
        C0: Optional[Union[float, np.ndarray]] = None,
        u: Optional[Union[float, np.ndarray]] = 0,
        beta: Optional[float] = None,
        mass: Optional[float] = 1.0) -> np.ndarray:
    r"""Calculate the covariance matrix :math:`\langle \mathbf{v}(0) \mathbf{v}(t)^{\intercal} \rangle` where :math:`\mathbf{v}(t)^{\intercal}  = (q, p, \mathbf{s}^{\intercal})` are the position and momentum of a harmonic oscillator + the auxiliary variable coordinates. See https://doi.org/10.1021/ct900563s for the notation.

    :param t: time(s) for which to compute the covariance matrix
    :type t: numpy.ndarray
    :param omega: frequency of the harmonic oscillator
    :type omega: float
    :param Ap: the momentum + auxiliary variable block of the drift matrix
    :type Ap: numpy.ndarray
    :param Cp: this parametrises the random forces; if omitted, the equilibrium thermal distribution is imposed
    :type Cp: numpy.ndarray, optional
    :param C0: Covariance matrix at time 0 (see notes)
    :type C0: numpy.ndarray, optional
    :param u: a second time variable, such that :math:`\langle \mathbf{v}(s) \mathbf{v}(t)^{\intercal} \rangle` is returned
    :type u: numpy.ndarray, optional
    :param beta: reciprocal temperature :math:`1/k_B T`
    :type beta: float, optional
    :param mass: mass of the harmonic oscillator
    :type mass: float, optional
    :return: covariance matrix
    :rtype: numpy.ndarray

    .. note::
        * If `C0` is not specified, the equilibrium distribution is assumed at :math:`t = 0`.
        * If a float is supplied, this is assumed to be the position variance, mass-weighted momenta are assumed to have a diagonal covariance :math:`1/\beta`
        * If a 2x2 array is supplied, this is assumed to be the position-momentum covariance, auxilliary momenta are assumed to have a diagonal covariance :math:`1/\beta`
        * Otherwise, the code expects a full (n+1,n+1) array, where `n` is the number of rows/columns in `Ap`.
    """
    
    from scipy.linalg import expm
    t = np.asarray(t)
    if t.ndim == 1:
        t = np.reshape(t, (1,-1,1,1)) 
    elif t.ndim == 2:
        t = t[...,None,None]
    else:
        raise RuntimeError('t should either be a one- or two-dimensional array')
    if np.any(t < 0):
        raise RuntimeError("Can only request positive times!")
    u = np.atleast_1d(u)
    if u.ndim != 1:
        raise RuntimeError('s should be a scalar or a one-dimensional array')
    u = np.reshape(u, (-1,1,1,1))
    if np.any(u < 0):
        raise RuntimeError("Can only request positive times!")
    if (beta is None and Cp is None) or (beta is not None and Cp is not None):
        raise RuntimeError(f"""Need to specify either beta (for sampling with detailed balance) or C_p (for a more general thermostat)""")
    Ap = np.atleast_2d(Ap)
    A_qp = get_Aqp(Ap, omega)
    n = Ap.shape[0]
    if Cp is None: # detailed balance fulfilled
        C_qp = get_Cqp_canonical(beta, n, omega)
    else: # frequency-dependent thermostat, e.g. quantum thermostat
        Cp = np.atleast_2d(Cp)
        # the first element of C_p is optimized to be kT for the quantum thermostat 
        # (need temperature in case of non-equilibrium starting covariance matrix)
        beta = 1/Cp[0, 0]
        D_qp = get_Dqp(Ap, Cp)
        C_qp = get_Cqp(A_qp, D_qp)
    if C0 is None: # at time zero have steady state distribution 
        term0 = 0
    else:
        C0_ = np.atleast_2d(C0)
        check_matrix(C0_, 'C_qp0')
        nqp = C0_.shape[0]
        C0 = np.eye(n+1) * (1/beta)
        if nqp in {1,2}:
            C0[:nqp,:nqp] = C0_[:nqp,:nqp]
        else:
            C0[:,:] = C0_
        term0 = expm(-A_qp * t) @ (C0 - C_qp) @ expm(-np.transpose(A_qp) * u)
    term1 = expm(-A_qp * np.abs(t - u)) @ C_qp
    tcf = term0 + term1
    sqm = np.sqrt(mass) # rescaling for mass != 1 
    sqmvec = np.ones(n+1)
    sqmvec[0] = 1/sqm
    sqmvec[1] = sqm
    sqmmat = sqmvec[:,None] * sqmvec[None,:]
    return np.squeeze(tcf * sqmmat)

def check_matrix(
        M : Union[float, np.ndarray],
        name : str) -> None:
    """
    Check if matrix M is 2-dimensional and square. If it fails a check raise a
    RuntimeError.

    Args:
        M : any if the input matrices
        name : the name of the matrix
    """
    if M.ndim != 2:
        raise RuntimeError(f"2-dimensional array expected for " + name + f", instead {M.ndim = }")
    if M.shape[0] != M.shape[1]:
        raise RuntimeError(f"Square " + name + f" expected, instead got {M.shape = }")  
    return

def get_Aqp(
        Ap : Union[float, np.ndarray], 
        omega : float) -> np.ndarray:
    """
    :param Ap: momentum + auxiliary variable block of the drift matrix, `shape(n, n)`
    :type Ap: numpy.ndarray
    :param omega: harmonic oscillator frequency 
    :type omega: float
    :return: full drift matrix, `shape(n+1, n+1)`
    :rtype: numpy.ndarray
    """
    check_matrix(Ap, 'Ap') 
    n = Ap.shape[0]
    Aqp = np.zeros((n+1, n+1))
    Aqp[1:,1:] = Ap
    Aqp[0, 1] = -1
    Aqp[1, 0] = omega**2
    return Aqp

def get_Cqp_canonical(
        beta : float, 
        n : int, 
        omega : float) -> np.ndarray: 
    r"""
    Compute the cacnonical covariance matrix for mass-weighted phase-space variables

    :param beta: reciprocal temperature :math:`1/k_B T`
    :type beta: float
    :param n: size of the `Ap` block of the drift matrix (i.e. 1 + number of auxiliary variables)
    :type n: int
    :param omega: harmonic oscillator frequency
    :type omega: float
    :return: covariance matrix for the canonical ensemble, `shape(n+1, n+1)`
    :rtype: numpy.ndarray
    """
    Cqp = np.eye(n+1) / beta # variance for momenta
    Cqp[0,0] /= omega**2 # variance for position
    return Cqp   

def get_Dqp(
        Ap : Union[float, np.ndarray], 
        Cp : Union[float, np.ndarray]) -> np.ndarray:
    """
    :param Ap: the momentum + auxiliary variable block of the drift matrix, `shape(n, n)`.
    :type Ap: numpy.ndarray
    :param Cp: the corresponding :math:`\mathbf{C}_p` matrix
    :type Cp: numpy.ndarray
    :return: :math:`\mathbf{D}_{qp}`-matrix [Eq (12) of https://doi.org/10.1021/ct900563s], `shape(n+1, n+1)`.
    :rtype: numpy.ndarray
    """
    check_matrix(Cp, 'C_p')
    n = Cp.shape[0]
    Dp = Ap @ Cp + Cp @ np.transpose(Ap)
    Dqp = np.zeros((n+1, n+1))
    Dqp[1:, 1:] = Dp
    return Dqp   


def get_Cqp(
        Aqp : Union[float, np.ndarray], 
        Dqp : Union[float, np.ndarray]) -> np.ndarray:
    """
    Compute the stationary covariance matrix for mass-weighted variables that solves

    .. math::
        \mathbf{A}_{qp} \mathbf{C}_{qp} + \mathbf{C}_{qp} \mathbf{A}_{qp}^{\intercal} = \mathbf{B}_{qp} \mathbf{B}_{qp}^{\intercal} = \mathbf{B}_{qp}^{\intercal}

    (see Appendix D of https://doi.org/20.500.11850/152344 for the algorithm)

    :param Aqp: full drift matrix, `shape(n+1, n+1)`
    :type Aqp: numpy.ndarray
    :param Dqp: :math:`\mathbf{B}_{qp} \mathbf{B}_{qp}^{\intercal}` where :math:`\mathbf{B}_{qp}` is the diffusion matrix, `shape(n+1, n+1)`
    :type Dqp: numpy.ndarray
    :return: full stationary covariance matrix, `shape(n+1, n+1)`
    :rtype: numpy.ndarray
    """
    import scipy.linalg as sclin
    eigs, O = sclin.eig(Aqp) # since A_qp is not symmetric eigs and O are complex
    O_inv = np.linalg.inv(O)
    n = Aqp.shape[0]
    Cqp = np.zeros_like(Aqp, dtype=complex)
    M = O_inv @ Dqp @ np.transpose(O_inv)
    for i in range(n): 
        # calculate each i, j element of C_qp individually
        for j in range(n):
            C = 0j
            for k in range(n):
                for l in range(n):
                    o = O[i, k] * M[k, l] * O[j, l]
                    a = eigs[k] + eigs[l]
                    C += o/a
            Cqp[i, j] = C
    return np.real(Cqp)