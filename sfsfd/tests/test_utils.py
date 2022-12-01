import numpy as np
from sfsfd.utils import fourier_to_polar, polar_to_fourier

def test_polar_conversion():
    """ Simple unit test for coeffs -> angle embedding and vice versa """

    # Embed and extract a diagonal vector
    x1 = np.ones(6) / np.sqrt(6.0)
    theta1 = fourier_to_polar(x1)
    assert(np.linalg.norm(polar_to_fourier(theta1) - x1) < 1.0e-8)
    # Extract and embed the first basis vector
    theta2 = np.ones(3) * 2.0 * np.pi
    x2 = polar_to_fourier(theta2)
    assert(np.linalg.norm(x2 - np.eye(4)[0]) < 1.0e-8)
    assert(np.linalg.norm(theta2 - fourier_to_polar(x2)) < 1.0e-8 or
           np.linalg.norm(fourier_to_polar(x2)) < 1.0e-8)
    # Extract and embed the first basis vector again
    theta3 = np.zeros(11)
    x3 = polar_to_fourier(theta3)
    assert(np.linalg.norm(x3 - np.eye(12)[0]) < 1.0e-8)
    assert(np.linalg.norm(theta3 - fourier_to_polar(x3)) < 1.0e-8 or
           np.linalg.norm(np.ones(11) * 2.0 * np.pi - fourier_to_polar(x2))
           < 1.0e-8)


if __name__ == "__main__":
    " Run the unit test manually when called as main "

    test_polar_conversion()
