import numpy as np

def radial_estimator(delta, qP_q1, sigP):
    """
    Calculate the radial estimator for the end-to-end distance

    Args:
    delta (float one dimensional array): distances for which to 
    evaluate the estimator
    qP_q1 (float or one dimensional array): sampled distances
    of the last and first bead
    sigP (float): standard deviation of the Gaussians, 
    sigP = sqrt(hbar^2 beta / (P m))

    Returns:
    N (one dimensional array): radial estimator as a 
    function of delta (same length), is averaged over all qP_q1
    distances
    """
    D, q = np.meshgrid(np.asarray(delta), np.asarray(qP_q1))
    #delta in rows and qP_q1 in columns
    a = (2 * np.pi * sigP**2)**(-0.5) / (D * q)
    e1 = np.exp(-((D - q)**2)/(2 * sigP**2))
    e2 = np.exp(-((D + q)**2)/(2 * sigP**2))
    X = a * (e1 - e2)
    N = np.mean(X, axis=0) # average over qP_q1 values in columns
    return N

def integral_del(x):
    xi = np.asarray(x)
    if xi.ndim==1:
        a = None
        return (np.roll(xi, -1, axis=a) - xi)[:-1]
    elif xi.ndim==2:
        a = 0
        return (np.roll(xi, -1, axis=a) - xi)[:-1, :-1]
    else:
        raise Exception('Wrong dimension {d} for x'.format(d=str(xi.ndim)))

def integral_av(y):
    yi = np.asarray(y)
    if yi.ndim==1:
        a = None
        return (np.roll(yi, -1, axis=a) + yi)[:-1]/2
    elif yi.ndim==2:
        a = 1
        return (np.roll(yi, -1, axis=a) + yi)[:-1, :-1]/2
    else:
        raise Exception('Wrong dimenstion {d} for x'.format(d=str(yi.ndim)))

def integrate(dx, y):
    """
    Integrate function y

    Args:
    dx (float, 1d or 2d array): distances between integration
    variable, if float -> equally spaced grid points, else
    the shape has to be the same as for y
    y (1d or 2d array): function to integrate, 1d for calculation of 
    one integral, 2d if array of integrals should be calculated

    for 2d arrays: obtain x and y as meshgrid, integration is done
    over axis 0
    """
    yi = np.asarray(y)
    try:
        l = len(dx)
    except:
        l = 0
    if l > 0 and dx.shape != y.shape:
        raise Exception('dx and y need to have the same shape! dx shape '
        'is {s_dx} while y shape is {s_y}'.format(s_dx = str(dx.shape), s_y=str(y.shape)))
    if y.ndim == 1:
        I = np.sum(dx*yi)
    elif y.ndim == 2:
        I = np.sum(dx*yi, axis=0)
    else:
        raise Exception('Wrong dimenstion {d} for y'.format(d=str(y.ndim)))
    return I


def radial_momentum(delta, N, p, hbar=1):
    """
    Calculate the spherically averaged momentum distribution

    Args:
    delta (float or one dimensional array): distances for which 
    the end-to-end distance is evaluated
    N (one dimensional array): radial end-to-end estimator, same
    length as delta
    p (float or one dimensional array): momenta for which to calculate
    the distribution
    """
    P, D = np.meshgrid(np.asarray(p), np.asarray(delta))
    P, Nm = np.meshgrid(np.asarray(p), np.asarray(N))
    # p in rows and delta dependence in columns
    integrand = Nm * D**2 * hbar / (P * D) * np.sin(P * D / hbar)
    del_D = integral_del(D)
    av_int = integral_av(integrand)
    I = integrate(del_D, av_int)
    p_new = integral_av(p)
    n = 4 * np.pi * I
    return p_new, n

def classical_momentum(beta, m, p):
    a = 4 * np.pi * (beta/(2 * np.pi *m))**(3/2)
    n = a * np.exp(- beta * p**2 /(2 * m))
    return n

def get_momentum(p, i_atom):
    """
    get momentum distribution from sampled momentum

    Args:
    p (3 dimensional numpy array): momentum of shape (Nsamples, natoms, 3)
    i_atom (list): indices of atoms for which to calculate the momentum distribution
    (to include all atoms of the same type)
    """
    Nsamples = p.shape[0]
    norm_p = np.zeros(Nsamples*len(i_atom))
    for j, i in enumerate(i_atom):
        norm_p[j*Nsamples:(j+1)*Nsamples] = np.linalg.norm(p[:, i], axis=1)
    return norm_p

