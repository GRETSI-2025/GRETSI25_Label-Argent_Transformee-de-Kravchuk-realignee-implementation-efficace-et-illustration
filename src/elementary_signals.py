# LOAD NECESSARY PYTHON LIBRARIES

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import cmocean


# HANDLE APPEARANCE OF THE PLOTS

my_cmap = cmocean.cm.cmap_d["deep"]
my_blue = my_cmap(0.6)
my_blue = my_blue[:3]
my_green = my_cmap(0.2)
my_green = my_green[:3]


# DETERMINISTIC SIGNALS


def the_chirp(N, observe=20, duration=15):
    """
    Chirp of unit energy for the 2-norm.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - observe (float): total duration of the observation window (in seconds).
        - duration (float): duration of the phenomenon of interest (in seconds).

    Returns:
        - signal (numpy.ndarray): discrete signal corresponding to the regular sampling of a chirp.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    time_t = np.linspace(-observe, observe, N + 1)
    envelop = np.zeros(N + 1)
    envelop[np.abs(time_t) < duration] = np.exp(
        -25 / (duration**2 - time_t[np.abs(time_t) < duration] ** 2)
    )
    freq = 0.5 + (time_t + duration) * (1 - 0.5) / (2 * duration)
    chirp = np.sin(2 * np.pi * freq * time_t)
    signal = chirp * envelop
    signal /= np.linalg.norm(signal, ord=2)

    return signal, time_t


def the_sine(N, observe=20, frequency=1):
    """
    Pure sine of unit energy for the 2-norm.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - observe (float): total duration of the observation window (in seconds).
        - frequency (float): frequency of the pure sine (in Hertz).

    Returns:
        - signal (numpy.ndarray): discrete signal corresponding to the regular sampling of a pure sine.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    time_t = np.linspace(-observe, observe, N + 1)
    signal = np.sin(2 * np.pi * frequency * time_t)
    signal /= np.linalg.norm(signal, ord=2)

    return signal, time_t


def the_dirac(N, observe=20, location=0):
    """
    Dirac mass of unit energy for the 2-norm.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - observe (float): total duration of the observation window (in seconds).
        - location (float): location of the Dirac mass (in seconds).

    Returns:
        - signal (numpy.ndarray): discrete signal corresponding to a Dirac mass at a fixed position.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    time_t = np.linspace(-observe, observe, N + 1)
    signal = np.zeros(np.shape(time_t))
    loc = np.argwhere(np.abs(time_t - location) == np.min(np.abs(time_t - location)))
    signal[loc[0]] = 1
    signal /= np.linalg.norm(signal, ord=2)

    return signal, time_t


def the_coherent_state(N, observe=20, theta=np.pi / 2, phi=0):
    """
    Coherent state of unit energy for the 2-norm.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - observe (float): total duration of the observation window (in seconds).
        - theta (float): polar angle indicating the position of the coherent state on the sphere (between 0 and pi).
        - phi (float): aligned azimuthal angle indicating the position of the coherent state on the sphere (between -pi and pi).

    Returns:
        - signal (numpy.ndarray): discrete signal corresponding to the regular sampling of a coherent state.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    if theta < 0 or theta > np.pi:
        error("")

    time_t = np.linspace(-observe, observe, N + 1)
    ns = np.arange(N + 1)
    signal = (
        np.sqrt(N + 1)
        / (4 * np.pi)
        * np.exp(
            1 / 2 * special.loggamma(N + 1)
            - 1 / 2 * special.loggamma(ns + 1)
            - 1 / 2 * special.loggamma(N - ns + 1)
        )
        * np.cos(theta / 2) ** (N - ns)
        * np.sin(theta / 2) ** (ns)
        * np.exp(-1j * (phi + np.pi) * ns)
    )
    signal /= np.linalg.norm(signal, ord=2)

    return signal, time_t


# PURE NOISE


def the_white_noise(N):
    """
    Complex white Gaussian noise normalized for the 2-norm.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.

    Returns:
        - wnoise (numpy.ndarray): discrete white noise generated as a standard Gaussian vector of size N+1.

    """

    wnoise = (np.random.randn(N + 1) + 1j * np.random.randn(N + 1)) / np.sqrt(2)
    wnoise /= np.linalg.norm(wnoise, ord=2)

    return wnoise


# NOISY SIGNALS


def the_noisy_chirp(N, snr=1, observe=20, duration=15, disp=False):
    """
    Chirp of unit energy for the 2-norm with additive Gaussian noise of unit energy for the 2-norm with fixed signal-to-noise ratio.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - snr (float): signal-to-noisy ratio (nonnegative, if Inf no noise).
        - observe (float): total duration of the observation window (in seconds).
        - duration (float): duration of the phenomenon of interest (in seconds).
        - disp (boolean, optional): if true display the generated signal.

    Returns:
        - nsignal (numpy.ndarray): discrete signal corresponding to the regular sampling of a noisy chirp.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """
    signal, time_t = the_chirp(N, observe, duration)
    wnoise = the_white_noise(N)

    if snr > 0:
        nsignal = signal + 1 / snr * wnoise
    else:
        nsignal = wnoise

    if disp:
        display_signal(nsignal, time_t)

    return nsignal, time_t


def the_noisy_sine(N, snr=1, observe=20, frequency=1, disp=False):
    """
    Pure sine of unit energy for the 2-norm with additive Gaussian noise of unit energy for the 2-norm with fixed signal-to-noise ratio.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - snr (float): signal-to-noisy ratio (nonnegative, if Inf no noise).
        - observe (float): total duration of the observation window (in seconds).
        - frequency (float): frequency of the pure sine (in Hertz).
        - disp (boolean, optional): if true display the generated signal.

    Returns:
        - nsignal (numpy.ndarray): discrete signal corresponding to the regular sampling of a noisy sine.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    signal, time_t = the_sine(N, observe, frequency)

    wnoise = the_white_noise(N)

    if snr > 0:
        nsignal = signal + 1 / snr * wnoise
    else:
        nsignal = wnoise

    if disp:
        display_signal(nsignal, time_t)

    return nsignal, time_t


def the_noisy_dirac(N, snr=1, observe=20, location=0, disp=False):
    """
    Dirac mass of unit energy for the 2-norm with additive Gaussian noise of unit energy for the 2-norm with fixed signal-to-noise ratio.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - snr (float): signal-to-noisy ratio (nonnegative, if Inf no noise).
        - observe (float): total duration of the observation window (in seconds).
        - location (float): location of the Dirac mass (in seconds).
        - disp (boolean, optional): if true display the generated signal.

    Returns:
        - nsignal (numpy.ndarray): discrete signal corresponding to a Dirac mass at a fixed position with additive Gaussian noise.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    signal, time_t = the_dirac(N, observe, location)

    wnoise = the_white_noise(N)

    if snr > 0:
        nsignal = signal + 1 / snr * wnoise
    else:
        nsignal = wnoise

    if disp:
        display_signal(nsignal, time_t)

    return nsignal, time_t


def the_noisy_coherent_state(N, snr=1, observe=20, theta=np.pi / 2, phi=0, disp=False):
    """
    Coherent state of unit energy for the 2-norm with additive Gaussian noise of unit energy for the 2-norm with fixed signal-to-noise ratio.

    Args:
        - N (integer): the number of samples in the generated signal will be N+1.
        - snr (float): signal-to-noisy ratio (nonnegative, if Inf no noise).
        - observe (float): total duration of the observation window (in seconds).
        - theta (float): polar angle indicating the position of the coherent state on the sphere (between 0 and pi).
        - phi (float): aligned azimuthal angle indicating the position of the coherent state on the sphere (between -pi and pi).
        - disp (boolean, optional): if true display the generated signal.

    Returns:
        - nsignal (numpy.ndarray): discrete signal corresponding to the regular sampling of a noisy coherent state.
        - time_t (numpy.ndarray): vector of time stamps at which the signal is sampled.

    """

    signal, time_t = the_coherent_state(N, observe, theta, phi)

    wnoise = the_white_noise(N)

    if snr > 0:
        nsignal = signal + 1 / snr * wnoise
    else:
        nsignal = wnoise

    if disp:
        display_signal(nsignal, time_t)

    return nsignal, time_t


# DISPLAY THE SIGNALS


def display_signal(nsignal, time_t=np.array([]), yticks=True):
    """
    Display the real part of a signal with respect to time.

    Args:
        - nsignal (numpy.ndarray): discrete signal, noisy or not, possibly complex valued.
        - time_t (numpy.ndarray, optional): vector of time stamps at which the signal is sampled.
        - yticks (boolean, optional): if true display the tick labels on y-axis.

    Returns:
        - time_t (numpy.ndarray): vector of time stamps at which the signal is displayed.

    """

    if len(time_t) == 0:
        time_t = np.arange(len(nsignal))

    plt.subplots(figsize=(5, 3))
    plt.plot(time_t, nsignal.real, color=my_blue, linewidth=1)
    plt.grid()
    plt.xlabel("$t$ (s)")
    if yticks == False:
        a = plt.gca()
        ya = a.axes.get_yticks()
        plt.yticks(ya[1:-1], [])

    plt.tight_layout()

    return time_t
