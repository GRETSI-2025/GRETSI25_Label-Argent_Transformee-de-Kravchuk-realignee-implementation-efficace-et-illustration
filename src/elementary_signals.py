# LOAD NECESSARY PYTHON LIBRARIES

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import cmocean

# LOAD FUNCTIONS TO PERFORM EXPERIMENTS ON ELEMENTARY SIGNALS

import time
from pkravchuk_transform import the_ptransform
from ssht_transform import the_new_transform, the_inverse_transform, the_spherical_angles

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
        raise NameError("The polar angle theta should be between 0 and pi.")

    if phi < -np.pi or theta > np.pi:
        raise NameError("The azimuthal angle phi should be between -pi and pi.")

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


# EXPERIMENTS TO EVALUATE COMPUTING TIME AND PRECISION


def time_and_precision(N = [256], R = 1, observe = 20, type = 'chirp', param = 15, snr = np.inf,calc = 'time-and-precision'):
    """
    Evaluate the computing time and reconstruction precision of the original and aligned transforms depending on signal size.

    Args:
        - N (list of integers, optional): sizes of signals to be explored.
        - R (integer,optional): number of realizations.
        - observe (float, optional): total duration of the observation window (default: 20 seconds).
        - type (string, optional): type of signal to be considered (default: chirp).
        - param (float, optional): parameter of the signal (default: duration of 15 seconds).
        - snr (float, optional): signal-to-noisy ratio (nonnegative, if Inf no noise).
        - calc (string, optional): either compute only time, only precision or both time and precision

    Returns:
        - to_signal (numpy.ndarray): mean and confidence interval of computing time of the original transform
        - tn_signal (numpy.ndarray): mean and confidence interval of computing time of the new aligned transform
        - ser_signal (numpy.ndarray): mean and confidence interval of signal-to-error ratio if calc indicates to compute time and precision

    """

    if calc == 'time':
        dotime  =  True
        doprec  =  False
    elif calc == 'precision':
        dotime  =  False
        doprec  =  True
    elif calc == 'time-and-precision':
        dotime  =  True
        doprec  =  True
    else:
        raise NameError("Type of signal not implemented")

        
    # computational time
    torigin = np.zeros((len(N),R))
    tnew    = np.zeros((len(N),R))

    if doprec :
        # precision
        prec    = np.zeros((len(N),R))

    for (i,n) in enumerate(N):

        # common grid of spherical angles on which to compute the Kravchuk transform
        (thetas, phis) = the_spherical_angles(n)

        if type == 'Dirac':
                
            for r in range(R):

                
                # generate signal
                signal,_     = the_noisy_dirac(n,observe = observe, location = param, snr = snr)

                if dotime:
                    if n <= 1024:
                        # compute the original and aligned Kravchul transform
                        ti_origin    = time.time()  
                        Ks, _,_      = the_ptransform(signal,1/2,thetas[:-1],phis+np.pi)
                        tf_origin    = time.time()  
        
                        # store computational time
                        dt_origin    = tf_origin - ti_origin
                        torigin[i,r] = dt_origin

                # compute the novel aligned Kravchuk transform of the signal
                ti_new       = time.time()  
                Fs,_,_       = the_new_transform(signal)
                tf_new       = time.time()  
                dt_new       = tf_new - ti_new
                tnew[i,r]    = dt_new

                if doprec :
                    # compute the inverse of the Kravchuk transform
                    isignal      = the_inverse_transform(Fs,n)
        
                    # compute the reconstruction error
                    prec[i,r]    = 1 / np.linalg.norm(isignal-signal,ord = 2)

            if R > 1:
                print('Finished R = ' + str(R) + ' realizations of signals of size N = ' + str(n) + '.')
            else:
                print('Finished R = ' + str(R) + ' realization of signals of size N = ' + str(n) + '.')
        
        elif type == 'sine':

            for r in range(R):
                
                # generate signal
                signal,_     = the_noisy_sine(n,observe = observe, frequency = param, snr = snr)

                if dotime:
                    if n <= 1024:
                        # compute the original and aligned Kravchul transform
                        ti_origin    = time.time()  
                        Ks, _,_      = the_ptransform(signal,1/2,thetas[:-1],phis+np.pi)
                        tf_origin    = time.time()  
        
                        # store computational time
                        dt_origin    = tf_origin - ti_origin
                        torigin[i,r] = dt_origin

                # compute the novel aligned Kravchuk transform of the signal
                ti_new       = time.time()  
                Fs,_,_       = the_new_transform(signal)
                tf_new       = time.time()  
                dt_new       = tf_new - ti_new
                tnew[i,r]    = dt_new

                if doprec :
                    # compute the inverse of the Kravchuk transform
                    isignal      = the_inverse_transform(Fs,n)
    
                    # compute the reconstruction error
                    prec[i,r]    = 1 / np.linalg.norm(isignal-signal,ord = 2)

            if R > 1:
                print('Finished R = ' + str(R) + ' realizations of signals of size N = ' + str(n) + '.')
            else:
                print('Finished R = ' + str(R) + ' realization of signals of size N = ' + str(n) + '.')
            
        elif type == 'chirp':

            for r in range(R):
                
                # generate signal
                signal,_     = the_noisy_chirp(n,observe = observe, duration = param, snr = snr)

                if dotime:
                    if n <= 1024:
                        # compute the original and aligned Kravchul transform
                        ti_origin    = time.time()  
                        Ks, _,_      = the_ptransform(signal,1/2,thetas[:-1],phis+np.pi)
                        tf_origin    = time.time()  
        
                        # store computational time
                        dt_origin    = tf_origin - ti_origin
                        torigin[i,r] = dt_origin

                # store computational time
                dt_origin    = tf_origin - ti_origin
                torigin[i,r] = dt_origin

                # compute the novel aligned Kravchuk transform of the signal
                ti_new       = time.time()  
                Fs,_,_       = the_new_transform(signal)
                tf_new       = time.time()  
                dt_new       = tf_new - ti_new
                tnew[i,r]    = dt_new

                if doprec :
                    # compute the inverse of the Kravchuk transform
                    isignal      = the_inverse_transform(Fs,n)
    
                    # compute the reconstruction error
                    prec[i,r]    = 1 / np.linalg.norm(isignal-signal,ord = 2)

            if R > 1:
                print('Finished R = ' + str(R) + ' realizations of signals of size N = ' + str(n) + '.')
            else:
                print('Finished R = ' + str(R) + ' realization of signals of size N = ' + str(n) + '.')
        
        else:
        
            raise NameError("Type of signal not implemented")

    # store means and confidence intervals
    if dotime:
        to_signal  = np.zeros((2,len(N)))
        tn_signal  = np.zeros((2,len(N)))
    if doprec :
        ser_signal = np.zeros((2,len(N)))
    
    # compute the mean over all realizations
    if dotime:
        to_signal[0,:]  = np.mean(torigin,axis = 1)
        tn_signal[0,:]  = np.mean(tnew,axis = 1)
    if doprec :
        ser_signal[0,:] = np.mean(prec,axis = 1)

    if R > 1:
        if dotime:
            # compute the 95% Gaussian confidence interval over all realizations
            to_signal[1,:]  = 1.96 / np.sqrt(R) *np.std(torigin,axis = 1)
            tn_signal[1,:]  = 1.96 / np.sqrt(R) *np.std(tnew,axis = 1)
        if doprec :
            ser_signal[1,:] = 1.96 / np.sqrt(R) *np.std(prec,axis = 1)

    if calc == 'time-and-precision':
        return to_signal, tn_signal, ser_signal
    elif calc == 'precision':
        return ser_signal
    elif calc == 'time':
        return to_signal, tn_signal