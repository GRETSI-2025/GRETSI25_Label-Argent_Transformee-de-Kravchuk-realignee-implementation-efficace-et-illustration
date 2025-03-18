import numpy                as np
import numpy.random         as npr
import matplotlib.pyplot             as plt
import os
from   scipy.io             import savemat, loadmat

from chirp_signals          import display_signal
from kravchuk_transform     import the_transform, the_zeros
from ssht_transform         import the_ssht_transform
from spherical_statistics   import the_distance, the_F_statistics, the_K_statistics, empirical_F, empirical_K
from detection_test         import the_test_statistic

def analyze_signal(signal, win = 1025, s = 32, folder = 'data/signals', name = 'chirp'):
    
    """
    analyze a long signal using overlapping windows of size win with a sliding strategy using a stride s
    """
    
    #signal: long signal to be analyzed
    #win:    size of the analysis window
    #s:      stride for the sliding window strategy
    #folder: name of the folder in which the Kravchuk spectrogram and spatial statistics are stored
    #name:   key for the signal used as prefix in the name of the file
    
    ## PREPARE STORAGE
    
    # create a folder to save the white noise samples
    if os.path.isdir('../'+folder) == False:
        os.makedirs('../'+folder)

    
    # PREPARE THE LOOP
    
    Nwin  = 0
    start = 0
    
    while start + win < len(signal):

        ## EXTRACT A WINDOW IN THE SIGNAL

        wsignal = signal[start:start + win]
    
        
        ## KRAVCHUK SPECTROGRAM ANALYSIS
        
        # compute the transform from Spin Spherical Harmonics
        Kb, thetas, phis   = the_ssht_transform(wsignal)
        
        # find the zeros
        zt, zp            = the_zeros(Kb,win-1,thetas,phis)

        
        ## SPATIAL STATISTICS ON ZEROS
        
        # explored range of r
        rs       = np.linspace(0,2*np.pi/np.sqrt(win-1),10**4)
        
        # compute the F-function
        Fb,_     = the_F_statistics(zt,zp,rs)
        
        # compute Ripley's K function
        Rb,_     = the_K_statistics(zt,zp)
        
        
        ## SAVE THE SPECTROGRAM AND THE SPATIAL FUNCTIONALS
        
        mdict    = {"win" : win, "s" : s, "wsignal" : wsignal,  "Kb" : Kb, "zt" : zt, "zp" : zp, "rs" : rs, "Fb" : Fb, "Rb" : Rb}

        savemat('../'+folder+'/'+name+'_'+str(Nwin)+'.mat',mdict)

        print('Analysis of window '+str(Nwin)+' done.')

        ## TRANSLATE TO NEXT WINDOW
        start    += s
        Nwin     += 1
        
    return



def extract_results(signal, time_t, win = 1025, s = 32, folder = 'data/signals', name = 'chirp', nfolder = 'data/noises', alpha = 0.05, m = 199, disp = False):
    
    """
    display the test performed on a long signal using overlapping windows of size win with a sliding strategy using a stride s
    """
    
    #signal: long signal to be analyzed
    #win:    size of the analysis window
    #s:      stride for the sliding window strategy
    #folder: name of the folder in which the Kravchuk spectrogram and spatial statistics are stored
    #name:   key for the signal used as prefix in the name of the file
    #alpha: confidence level of the test
    #m: number of samples under the null hypothesis to be used

    # DISPLAY THE SIGNAL AND A ZOOM ON ITS CENTRAL PART

    if disp:
        time_t     = display_signal(signal,time_t) 
    idc            = np.argwhere(np.abs(time_t)<0.05) 
    if disp:
        _          = display_signal(signal[idc],time_t[idc]) 

    
    ## SETTINGS OF THE TEST AND PREPARE STORAGE
    
    k              = int(alpha*(m+1))
    functional     = 'F'
    nr             = 10**4
    S_m            = np.zeros((m+1,nr))

    
    #LOAD THE FUNCTIONAL STATISTICS OF THE NOISE SAMPLES
    
    for n in np.arange(m):
        ndict         = loadmat(nfolder+'/noise_'+str(n)+'.mat',squeeze_me=True);
        if functional == 'F':
            S_m[n,:]  = ndict["Fn"]
        else:
            S_m[n,:]  = ndict["Rn"]


    # PREPARE THE LOOP
    
    Nwin  = 0
    start = 0
    T_exp   = []
    T_mk    = []

    while start + win < len(signal):

        ## LOAD THE SPECTROGRAM AND THE SPATIAL FUNCTIONALS
        
        mdict             = loadmat(folder+'/'+name+'_'+str(Nwin)+'.mat',squeeze_me=True)
        
                
        ## STORE THE SPATIAL FUNCTIONAL OF THE WINDOWED SIGNAL
        
        if functional     == 'F':
            S             = mdict["Fb"]
        else:
            S             = mdict["Rb"]
        S_m[-1,:]         = S


        ## COMPUTE AND STORE THE SUMMARY STATISTIC

        S0                = np.sum(S_m, axis = 0)/(m+1)
        T_m               = np.zeros(m)
        for i in np.arange(m):
                T_m[i]    = the_test_statistic(S_m[i,:],S0)
        T_m.sort()
        t_exp             = the_test_statistic(S,S0)
        T_exp.append(t_exp)
        T_mk.append(T_m[-k])
        

        ## TRANSLATE TO NEXT WINDOW

        start    += s
        Nwin     += 1
    

    ## DISPLAY THE RESULTS

    # location of the centers of the windows
    ctrs          = win//2 + s*np.arange(0,Nwin)
    time_T        = time_t[ctrs]

    if disp:    
        plt.plot(time_T,T_exp,'.k',label = 'summary statistic')
        plt.plot(time_T,T_mk,'-r',linewidth=2,label = 'detection threshold')
        plt.grid()
        plt.legend();
        plt.xlabel('time (sec)');
    
    
        detect = 100 * np.mean(np.array(T_exp)>=np.array(T_mk))
        print("Windows in which a signal is detected: %.2f%%.  Expected if pure noise: 5%%." % detect)
        
    return T_exp, time_T, T_mk






