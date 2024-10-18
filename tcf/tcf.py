import numpy as np
from typing import Optional

def autocorrelate(A: np.ndarray, axis: Optional[int] = 0) -> np.ndarray:
    """
    Compute the autocorrelation function of an array along the specified axis.
    
    :param A: array containing the time series
    :type A: np.ndarray
    :param axis: axis corresponding to time, defaults to 0
    :type axis: Optional[int], optional
    :return: autocorrelation function of the array input
    :rtype: np.ndarray
    """
    len_a = A.shape[axis]
    len_fft = 2 * len_a
    
    # Normalization factor (decreasing window size)
    norm_tcf = np.arange(len_a, 0, -1, dtype=int)
    
    # Expand normalization dimensions    
    norm_tcf = np.expand_dims(norm_tcf, tuple(i for i in range(A.ndim) if i != axis))
    
    # Compute FFT of A
    ftA = np.fft.rfft(A, axis=axis, n=len_fft)
    
    # Power spectrum (multiply by complex conjugate)
    ftA *= np.conj(ftA)
    
    # Inverse FFT to get the auto-correlation
    autocorr = np.fft.irfft(ftA, axis=axis, n=len_fft)
    
    # Slice the auto-correlation to original signal length
    autocorr = np.take(autocorr, np.arange(len_a), axis=axis)
    
    return autocorr / norm_tcf
