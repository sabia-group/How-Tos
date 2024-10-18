#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
Spectral densities, memory-friction kernels, harmonic discretisation schemes.
'''


from __future__ import print_function, division, absolute_import
import numpy as np
from typing import Union


class BaseSpectralDensity(object):

    def __init__(self,
                 mass: float,
                 Nmodes: int,
                 *args, 
                 **kwargs):
        """Model spectral densities for one-dimensional systems with bilinear coupling.

        :param mass: mass of the bath modes
        :type mass: float
        :param Nmodes: number of oscillators in the harmonic bath discretisation
        :type Nmodes: int
        """
        
        self.Nmodes = Nmodes
        self.c, self._frequencies, self.bath_mass = self.bath_params(mass)
        self.ww = self._frequencies**2
        self.mww = self.bath_mass*self.ww
    
    @property
    def frequencies(self):
        """Array of frequencies computed in the harmonic discretization of the spectral density.
        """
        return np.copy(self._frequencies)

    def quadpoints(self):
        raise NotImplementedError
    
    def bath_params(self, mass: float):
        m = np.ones(self.Nmodes)
        m *= mass
        w = self.quadpoints()
        kappa = np.sqrt(self.exact_reorganisation()/(2*self.Nmodes))
        c = kappa * np.sqrt(m) * w
        return c, w, m
    
    def quadrature(self, f):
        r"""Calculate the quadrature approximation to the integral

        .. math::
            \int_0^{\infty} J(\omega) f(\omega) \, \mathrm{d}\omega \approx 
            \frac{\pi}{2} \sum_{n=1}^{N_{\text{bath}}}  \frac{c_n^2}{m \omega_n} f(\omega_n)
          
        where :math:`J(\omega)` is the spectral density.

        :param f: input function evaluated at the quadrature points (`self.frequencies`).
        :f type: numpy.ndarray
        :return: integral over the input with the spectral density as the integration kernel
        :rtype: float
        """

        return (np.pi/2) * np.sum(self.c**2/(self.bath_mass*self._frequencies) * f, axis=-1)
    
    def l_quadrature(self, f):
        """Same as `quadrature` but using :math:`\Lambda(\omega) = J(\omega)/\omega` as the integration kernel.
        """

        return (np.pi/2) * np.sum(self.c**2/(self.bath_mass*self._frequencies**2) * f, axis=-1)
    
    def exact_reorganisation(self):
        raise NotImplementedError
    
    def reorganisation_energy(self):
        return (4/np.pi) * self.quadrature(1/self._frequencies)
    
    def J(self, omega: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Spectral density at frequency `omega`.
        """
        raise NotImplementedError
    
    def Lambda(self, omega: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Spectral density at frequency `omega`, divided by `omega`.
        """
        raise NotImplementedError
    
    def K(self, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Memory-friction kernel at time `t`.
        """
        raise NotImplementedError
    


class ExpOhmic(BaseSpectralDensity):

    def __init__(self, 
                 mass: float,
                 Nmodes: int,
                 eta: float, 
                 omega_cut: float, 
                 *args, **kwargs) -> None:
        """Exponentially damped Ohmic spectral density

        .. math::
            J(\omega) = \eta \omega \exp(-\omega/\omega_c)

        with discrete frequencies calculated according to Craig and Manolopoulos (2004), https://doi.org/10.1063/1.1850093.

        :param mass: mass of the bath modes
        :type mass: float
        :param Nmodes: number of oscillators in the harmonic bath discretisation
        :type Nmodes: int
        :param eta: static friction coefficient
        :type eta: float
        :param omega_cut: cut-off frequency
        :type omega_cut: float
        """
        self.eta = eta
        self.omega_cut = omega_cut
        super().__init__(mass, Nmodes, *args, **kwargs)

    def J(self, omega: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        return self.eta * np.abs(omega) * np.exp(-np.abs(omega)/self.omega_cut)
    
    def Lambda(self, omega: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        return self.eta * np.exp(-np.abs(omega)/self.omega_cut)
    
    def K(self, t: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        wc = self.omega_cut
        return 2 * self.eta * wc / (1 + (wc*t)**2) / np.pi
    
    def quadpoints(self) -> np.ndarray:
        return -self.omega_cut * np.log(
            (np.arange(1, self.Nmodes+1)-1/2) / self.Nmodes
        )[::-1]
    
    def exact_reorganisation(self) -> float:
        return (4/np.pi) * self.eta * self.omega_cut

