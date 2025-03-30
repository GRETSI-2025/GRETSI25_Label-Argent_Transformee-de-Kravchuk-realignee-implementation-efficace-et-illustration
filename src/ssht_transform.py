# LOAD NECESSARY PYTHON LIBRARIES

import numpy as np
import pyssht as ssht


def the_ssht_transform(x, method="MW", p=1 / 2):
    """
    Compute the p-Kravchuk transform of a discrete signal leveraging the connection with Spin Weighted Spherical Harmonics.

    Args:
        - x (numpy.ndarray): discrete signal, noisy or not, possibly complex valued.
        - method (string, optional): method of decomposition to use in pyssht (default "MW" for McEwen & Wiaux sampling).
        - p (float, optional): parameter of the transform between 0 and 1 (default 0.5 corresponding to the original transform).

    Returns:
        - f (numpy.ndarray): p-Kravchuk transform, complex-valued, evaluated at a discrete set of points on the sphere.
        - thetas (numpy.ndarray): polar angles at which the transform is computed.
        - phis (numpy.ndarray): azimuthal angles at which the transform is computed.
    """

    N = x.shape[0] - 1

    if np.mod(N, 2):
        print("Error: length of signal should be odd.")
    else:
        # largest index explored
        ell = int(N / 2)

        # band-limit
        L = ell + 1

        # spherical coordinates of the pixel centers used by the SSHT toolbox
        (thetas, phis) = ssht.sample_positions(L, Method=method)  # vectors
        (Thetas, Phis) = ssht.sample_positions(L, Method=method, Grid=True)  # meshgrid
        thetas = np.pi - thetas

        # store the spin-spherical harmonic coefficient
        flm = np.zeros(L**2, dtype="complex128")

        # compute the angle beta associated to p
        beta = 2 * np.arcsin(np.sqrt(p))

        for m in np.arange(-ell, ell + 1):
            n = ell + m
            index = ssht.elm2ind(ell, m)
            flm[index] = np.conj(x[n])

        slm = ssht.rotate_flms(flm, 0, beta, 0, L)
        f = ssht.inverse(slm, L, Spin=ell, Method=method)
        f *= np.sqrt(4 * np.pi / (N + 1))
        f *= np.exp(1j * N * Phis / 2)
        f *= (-1) ** (N / 2)

    return f, thetas, phis


def rotate_signal(x, theta=0, phi=0, method="MW"):
    """
    Rotate a signal leveraging the connection with Spin Weighted Spherical Harmonics.

    Args:
        - x (numpy.ndarray): discrete signal, noisy or not, possibly complex valued.
        - theta (float, optional): polar angle of the rotation of the signal (default 0).
        - phi (float, optional): azimuthal angle of the rotation of the signal (default 0).
        - method (string, optional): method of decomposition to use in pyssht (default "MW" for McEwen & Wiaux sampling).

    Returns:
        - y (numpy.ndarray): rotated signal.
    """

    if theta < 0 or theta > np.pi:
        raise NameError("The polar angle theta should be between 0 and pi.")

    if phi < -np.pi or theta > np.pi:
        raise NameError("The azimuthal angle phi should be between -pi and pi.")

    N = x.shape[0] - 1

    if np.mod(N, 2):
        print("Error: length of signal should be odd.")
    else:
        # largest index explored
        ell = int(N / 2)

        # band-limit
        L = ell + 1

        # store the spin-spherical harmonic coefficient
        flm = np.zeros(L**2, dtype="complex128")

        for m in np.arange(-ell, ell + 1):
            n = ell + m
            index = ssht.elm2ind(ell, m)
            flm[index] = np.conj(x[n])

        slm = ssht.rotate_flms(flm, phi, theta, 0, L)

        # store the reconstructed signal coefficients
        y = np.zeros((N + 1,), dtype="complex128")
        for m in np.arange(-ell, ell + 1):
            n = ell + m
            index = ssht.elm2ind(ell, m)
            y[n] = np.conj(slm[index])

    return y


def the_new_transform(x, method="MW"):
    """
    Compute the novel aligned Kravchuk transform of a discrete signal defined as an expansion on Spin Weighted Spherical Harmonics.

    Args:
        - x (numpy.ndarray): discrete signal, noisy or not, possibly complex valued.
        - method (string, optional): method of decomposition to use in pyssht (default "MW" for McEwen & Wiaux sampling).

    Returns:
        - f (numpy.ndarray): novel aligned Kravchuk transform, complex-valued, evaluated at a discrete set of points on the sphere.
        - thetas (numpy.ndarray): polar angles at which the transform is computed.
        - phis (numpy.ndarray): azimuthal angles at which the transform is computed.
    """

    N = x.shape[0] - 1

    if np.mod(N, 2):
        print("Error: N+1 should be even")
    else:
        # largest index explored
        ell = int(N / 2)

        # band-limit
        L = ell + 1

        # spherical coordinates of the pixel centers used by the SSHT toolbox
        (thetas, phis) = ssht.sample_positions(L, Method=method)  # vectors
        (Thetas, Phis) = ssht.sample_positions(L, Method=method, Grid=True)  # meshgrid

        # store the spin-spherical harmonic coefficients
        flm = np.zeros(L**2, dtype="complex128")
        for m in np.arange(-ell, ell + 1):
            n = ell + m
            index = ssht.elm2ind(ell, m)
            flm[index] = np.conj(x[n])

        # compute the inverse SSHT transform which coincides with the new Kravchuk transform
        f = ssht.inverse(flm, L, Spin=ell, Method=method)

        # center the frequencies
        phis += -np.pi

    return f, thetas, phis


def the_inverse_transform(f, N, method="MW"):
    """
    Compute the inverse of the novel aligned Kravchuk transform of a discrete signal defined as an expansion on Spin Weighted Spherical Harmonics.

    Args:
        - f (numpy.ndarray): novel aligned Kravchuk transform, complex-valued, evaluated at a discrete set of points on the sphere.
        - N (even integer): size of the original signal to be retrieved from its Kravchuk transform will be N+1.
        - method (string, optional): method of decomposition to use in pyssht (default "MW" for McEwen & Wiaux sampling).

    Returns:
        - x (numpy.ndarray): discrete signal having as aligned Kravchuk transform f.
    """

    if np.mod(N, 2):
        print("Error: N+1 should be even")
    else:
        # largest index explored
        ell = int(N / 2)

        # band-limit
        L = ell + 1

        # compute the SSHT transform which coincides with the new Kravchuk inverse transform
        flm = ssht.forward(f, L, Spin=ell, Method=method)

        # store the reconstructed signal coefficients
        x = np.zeros((N + 1,), dtype="complex128")
        for m in np.arange(-ell, ell + 1):
            n = ell + m
            index = ssht.elm2ind(ell, m)
            x[n] = np.conj(flm[index])

    return x


def the_spherical_angles(N, method="MW"):
    """
    Indicate the spherical angles at which the novel aligned Kravchuk transform of a signal is evaluated leveraging the pyssht library.

    Args:
        - N (even integer): size of the analyzed signal is N+1.
        - method (string, optional): method of decomposition to use in pyssht (default "MW" for McEwen & Wiaux sampling).

    Returns:
        - thetas (numpy.ndarray): polar angles at which the transform is computed.
        - phis (numpy.ndarray): azimuthal angles at which the transform is computed.
    """

    # largest index explored
    ell = int(N / 2)

    # band-limit
    L = ell + 1

    (thetas, phis) = ssht.sample_positions(L, Method=method)

    # adjust the azimuthal angle to fit the aligned transform definition
    phis = phis - np.pi

    return thetas, phis
