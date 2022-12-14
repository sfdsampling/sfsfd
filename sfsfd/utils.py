import numpy as np


def polar_to_fourier(theta):
    """ Convert polar angles into Fourier (wave function) coefficients.

    Args:
        theta (np.ndarray): A numpy array of length n - 1 containing polar
            angles.
    Returns:
        np.ndarray: A numpy array of length n containing the Fourier
        coefficients corresponding to the given polar angles.
    """

    # Initialize output array
    n = len(theta) + 1
    cx = np.ones(n)
    # Calculate the output array
    for i in range(n):
        for j in range(n - 1 - i):
            cx[i] *= np.cos(theta[j])
        if i > 0:
            cx[i] *= np.sin(theta[n - 1 - i])
    return cx


def fourier_to_polar(cx):
    """ Convert Fourier (wave function) coefficients into polar angles.
    Args:
        cx (np.ndarray): A numpy array of length n containing the Fourier
            coefficients corresponding to the given polar angles.

    Returns:
        theta (np.ndarray): A numpy array of length n - 1 containing polar
        angles.

    """

    # Initialize output array
    n = len(cx)
    theta = np.ones(n - 1)
    # Calculate the output array
    for i in range(n - 1):
        xx = cx[n - 1 - i]
        for j in range(i):
            xx /= np.cos(theta[j])
        theta[i] = np.arcsin(xx)
    return theta
