import numpy as np
import matplotlib.pyplot as plt
#import scipy.signal import chirp, find_peaks, peak_widths

def squarefunc(x):
    """
    Returns a squarefunction
    with value one if abs(x) < R,
    where R is a set value.
    """
    if np.abs(x) < R:
        return 1
    else:
        return 0

def W(k):
    """
    Function returning the fourier conjugate W(k)
    of squarefunc(x)
    """
    fac1 = np.sqrt(2 / np.pi)
    fac2 = np.sin(R*k) / k 
    return(fac1*fac2)

R = 0.5

k_lin = np.linspace(-100*R, 100*R, 2000)
W_hat = W(k_lin)

#Calculating the FWHM of the peak
#of W(k)

difference = np.max(W_hat) - np.min(W_hat)
half_max = difference / 2
nearest = np.argmin(np.abs(W_hat - half_max))
k_argmax = np.argmax(W_hat)
HWHM = np.abs(k_lin[k_argmax] - k_lin[nearest]) 
FWHM = 2*HWHM
print(FWHM)