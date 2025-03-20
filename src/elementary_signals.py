# LOAD NECESSARY PYTHON LIBRARIES

import numpy             as np
import matplotlib.pyplot as plt
import numpy.random      as npr
import cmocean 


# HANDLE APPEARANCE OF THE PLOTS

my_cmap  = cmocean.cm.cmap_d['deep']
my_blue  = my_cmap(0.6)
my_blue  = my_blue[:3]
my_green = my_cmap(0.2)
my_green = my_green[:3]


# DETERMINISTIC SIGNALS

def the_chirp(N,observe = 20,duration = 15):
    
    """
    chirp of unit energy for the 2-norm
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    
    time_t   = np.linspace(-observe, observe, N+1)
    envelop  = np.zeros(N+1)
    envelop[np.abs(time_t)<duration] = np.exp(-25/(duration**2-time_t[np.abs(time_t)<duration]**2))
    freq     = 0.5 + (time_t + duration)*(1 - 0.5)/(2*duration)
    chirp    = np.sin(2*np.pi*freq*time_t)
    signal   = chirp*envelop
    signal   /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t
    

def the_sine(N,observe = 20,frequency = 1):
    
    """
    chirp of unit energy for the 2-norm
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # frequency: frequency of the pure sine
    
    
    time_t   = np.linspace(-observe, observe, N+1)
    signal    = np.sin(2*np.pi*frequency*time_t)
    signal   /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t


def the_dirac(N,observe = 20,location = 0):
    
    """
    Dirac mass of unit energy localized at a given time
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # location: position of the Dirac mass
    
    
    time_t         = np.linspace(-observe, observe, N+1)
    signal         = np.zeros(np.shape(time_t))
    loc            = np.argwhere(np.abs(time_t - location) == np.min(np.abs(time_t - location)))
    signal[loc[0]] = 1
    signal        /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t
    
    

# PURE NOISE

def the_white_noise(N):
    
    """
    complex white Gaussian noise of length N+1 normalized for the 2-norm
    """
    
    # N+1: number of samples
    
    wnoise = (np.random.randn(N+1)+1j*np.random.randn(N+1))/np.sqrt(2)
    wnoise /= np.linalg.norm(wnoise,ord=2)
    
    return wnoise



# NOISY SIGNALS

def the_noisy_chirp(N,snr = 1,observe=20,duration=15,disp = False):
    
    """
    noisy chirp
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    signal, time_t = the_chirp(N,observe,duration)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t


def the_noisy_sine(N,snr = 1,observe=20,frequency=1,disp = False):
    
    """
    noisy sine
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # frequency: frequency of the pure sine
    
    signal, time_t = the_sine(N,observe,frequency)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t
    

def the_noisy_dirac(N,snr = 1,observe=20,location=0,disp = False):
    
    """
    noisy dirac
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # location: position of the Dirac mass
    
    signal, time_t = the_dirac(N,observe,location)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t
    

# DISPLAY THE SIGNALS

def display_signal(nsignal,time_t=np.array([]), yticks=True):
                   
    """
    display the real part of a signal along time
    """
    
    # nsignal: complex signal to be analyzed of length N+1
    # time_t: (optional) time range of observation (default: [0,1, ...,N])
    # yticks: (optional) if false remove the ticks on the y axis
    
    if len(time_t) == 0:
        time_t = np.arange(len(nsignal))

    plt.subplots(figsize=(5, 3))
    plt.plot(time_t,nsignal.real,color = my_blue, linewidth = 1);
    plt.grid()
    plt.xlabel('$t$ (s)')
    if yticks == False:
        a  = plt.gca()
        ya = a.axes.get_yticks()
        plt.yticks(ya[1:-1],[])
        
    plt.tight_layout()
    
    return time_t
    
