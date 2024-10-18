#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   special.py
@Time    :   2023/08/29 11:24:22
@Author  :   George Trenins
@Contact :   gstrenin@gmail.com
@Desc    :   Special functions
'''


from __future__ import print_function, division, absolute_import
import numpy as np
import scipy as sp


def sinc(x, nu=0):
    """Return sinc(x) = sin(x)/x or its nu'th derivative (nu = 0, 1, or 2)
    """
    y = x/np.pi
    if nu == 0:
        return sp.special.sinc(y)
    else:
        eps = 1.0e-4
        mask = np.abs(y) > eps
        temp = np.where(mask, x, 1.0)
        sinc_x = sp.special.sinc(y)
        if nu == 1:
            ans = np.where(mask, 
                           (np.cos(x)-sinc_x)/temp,
                           x * (x**2 * ((x**2/45360 - 1/840) * x**2 + 1/30) - 1/3)
                           )
        elif nu == 2:
            ans = np.where(mask, 
                           2*sinc_x/temp**2 - 2*np.cos(x)/temp**2 - sinc_x,
                           x**2 * ((x**2/6480 - 1/168) * x**2 + 1/10) - 1/3)
        else:
            raise NotImplementedError(f"Could not calculate the nu = {nu} derivative of sinc.")
        return ans


def logcosh(x):
    """Logarithm of cosh(x).
    """
    s = np.sign(x)
    x = s*x # real part now positive
    return x + np.log1p(np.exp(-2*x)) - np.log(2.0)

def logsinh(x):
    """Logarithm of sinh(x).
    """
    s = np.sign(x)
    x = s*x # real part now positive
    ans = np.where(np.real(s) < 0, 1j*np.pi, 0)
    ans += x + np.log1p(-np.exp(-2*x)) - np.log(2.0)
    return ans

    
if __name__ == "__main__":
    from gtlib.numderiv import gradient, hessian
    x = np.array([-2.0e-4, -3.1e-5, 0.0, 1.2e-6, 4.7e-3])
    print(sinc(x))
    print()
    print(sinc(x, nu=1))
    print([gradient(sinc, y, order=2, h=1.0e-6).item() for y in x])
    print()
    print(sinc(x, nu=2))
    print()
    print([gradient(lambda x: sinc(x, nu=1), y, order=2, h=1.0e-6).item() for y in x])
    n = 10
    re = np.random.uniform(1.0, 4, size=n)
    im = np.random.uniform(-4, 4, size=n)
    print("\n\n")
    for r in re:
        # for i in im:
        z = r #+ 1j*i
        print(np.exp(logcosh(z)) - np.cosh(z))
        print(np.exp(logsinh(z)) - np.sinh(z))
        print()