import numpy             as np
import matplotlib.pyplot as plt
import numpy.random      as npr
import cmocean 

sapin    = (0.0353, 0.3216, 0.1569)
my_cmap  = cmocean.cm.cmap_d['deep']
my_blue  = my_cmap(0.6)
my_blue  = my_blue[:3]
my_green = my_cmap(0.2)
my_green = my_green[:3]



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
    
    

def the_gaussian(N,observe = 20,mu = 0, sigma = 1):
    
    """
    Gaussian of unit energy for the 2-norm centered at mu and of width sigma
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # mu: location of the mode
    # sigma: width of the mode
    
    
    time_t   = np.linspace(-observe, observe, N+1)
    signal   = (np.pi*sigma**2)**(-1/4)*np.exp(-(time_t - mu)**2/(2*sigma**2))
    signal   /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t
    

def the_constant(N,observe = 20,value = 1):
    
    """
    Dirac mass of unit energy localized at a given time
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # value: value of the constant

    
    time_t         = np.linspace(-observe, observe, N+1)
    signal         = value * np.ones(np.shape(time_t))
    signal        /= np.linalg.norm(signal,ord=2)
        
    return signal, time_t


def the_white_noise(N):
    
    """
    complex white Gaussian noise of length N+1 normalized for the 2-norm
    """
    
    # N+1: number of samples
    
    wnoise = (np.random.randn(N+1)+1j*np.random.randn(N+1))/np.sqrt(2)
    wnoise /= np.linalg.norm(wnoise,ord=2)
    
    return wnoise
    

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
    

def the_noisy_gaussian(N,snr = 1,observe = 20,mu = 0, sigma = 1,disp = False):
    
    """
    noisy Gaussian
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # mu: location of the mode
    # sigma: width of the mode
    
    signal, time_t = the_gaussian(N,observe,mu,sigma)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t
    

def the_noisy_constant(N,snr = 1,observe = 20,value = 1,disp = False):
    
    """
    noisy constant
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # value: value of the constant
    
    signal, time_t = the_constant(N,observe,value)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t
    

def the_bird_call(N,observe=0.5,duration=0.05):
    
    """
    synthetic sample of a pure bird call
    """
    
    # N+1: number of samples
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    
    time_t   = np.linspace(-observe, observe, N+1)
    envelop  = np.zeros(N+1)
    envelop[np.abs(time_t)<duration] = np.exp(-0.005/(duration**2-time_t[np.abs(time_t)<duration]**2))
    freq     = 5000 + (time_t + duration)*(8000 - 5000)/(2*duration)
    chirp    = np.sin(2*np.pi*freq*time_t)
    signal   = chirp*envelop
    signal   /= np.linalg.norm(signal,ord=2)

        
    return signal, time_t


def the_noisy_bird(N,snr = 1,observe=0.5,duration=0.05,disp = False):
    
    """
    noisy bird call
    """
    
    # N+1: number of samples
    # snr: signal-to-noise ratio
    # observe: total duration of the measurement
    # duration: duration of the phenomenon of interest
    
    signal, time_t = the_bird_call(N,observe,duration)
    
    wnoise = the_white_noise(N)
    
    if snr > 0:
        nsignal = signal + 1/snr*wnoise
    else:
        nsignal = wnoise
        
    if disp:
        display_signal(nsignal,time_t)
        
    return nsignal, time_t


def display_signal(nsignal,time_t=np.array([]), yticks=True):
                   
    """
    display the real part of a signal along time
    """
    
    # nsignal: complex signal to be analyzed of length N+1
    # time_t: (optional) time range of observation (default: [0,1, ...,N])
    
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

def display_noisy_and_true(nsignal,signal,time_t=np.array([]),color = my_blue):
                   
    """
    display the real part of a signal along time
    """
    
    # nsignal: complex signal to be analyzed of length N+1
    # time_t: (optional) time range of observation (default: [0,1, ...,N])
    
    if len(time_t) == 0:
        time_t = np.arange(len(nsignal))

    plt.subplots(figsize=(5, 3))
    plt.plot(time_t,nsignal.real,color = 'gray', linewidth = 1,alpha = 0.5,label = r'observation $\boldsymbol{y}$')
    plt.plot(time_t,signal.real,color = my_blue, linewidth = 1.5,label = r'pure chirp $\boldsymbol{x}$')
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.grid()
    plt.xlabel('$t$ (s)')
    plt.xlim([time_t[0], time_t[-1]]);
    
    return time_t


def display_signal_and_zoom(nsignal,time_t,duration):
                   
    """
    display the real part of a signal along time and a zoom on the chirp
    """
    
    # nsignal: complex signal to be analyzed of length N+1
    # time_t: time range of observation
    # duration: length of the chirp

    # time stamps corresponding to the synthetic bird call
    idc            = np.argwhere(np.abs(time_t)<duration) 
    time_dc        = np.squeeze(time_t[idc])
    
    # bounds of the y-axis
    ymin = - 1.1 * np.max( np.abs(nsignal) ) 
    ymax = 1.1 * np.max( np.abs(nsignal) ) 
    
    
    fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout = True)
    ax1.fill_between(time_dc,ymin*np.ones(np.shape(time_dc)),ymax*np.ones(np.shape(time_dc)),color = my_green,alpha = 0.65)
    ax1.plot(time_t,nsignal.real,color = my_blue, linewidth = 1)
    ax1.grid()
    ax1.set_xlim([time_t[0], time_t[-1]]);
    ax1.set_ylim([ymin, ymax])
    ax2.plot(time_t[idc],nsignal.real[idc],color = my_blue, linewidth = 1)
    ax2.grid()
    ax2.set_xlim([time_t[idc[0]], time_t[idc[-1]]]);
    ax2.set_ylim([ymin, ymax])
    ax2.set_xlabel('$t$ (s)') ;
    
    
