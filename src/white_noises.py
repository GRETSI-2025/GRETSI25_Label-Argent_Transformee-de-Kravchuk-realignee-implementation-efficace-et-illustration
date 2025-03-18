import numpy                as np
import numpy.random         as npr
import os
from   scipy.io             import savemat, loadmat

from chirp_signals          import the_chirp, the_white_noise, the_noisy_chirp
from kravchuk_transform     import the_transform, the_zeros
from ssht_transform         import the_ssht_transform
from kravchuk_display       import signal_display, planar_display, spherical_display
from spherical_statistics   import the_distance, the_F_statistics, the_K_statistics, empirical_F, empirical_K

def noise_samples(N, m=199, folder = '../data/noises', start = 0):
    
    """
    generate m samples of complex white Gaussian noise and perform all required analysis for zero-based detection
    """
    
    #N+1: number of points in the white noise
    #m:   number of realizations for the Monte Carlo envelope test
    #folder: name of the folder in which the white noise samples are stored
    #start: index at which to start the numbering
    
    ## PREPARE STORAGE
    
    # create a folder to save the white noise samples
    if os.path.isdir(folder) == False:
        os.mkdir(folder)

    
    for n in np.arange(start,m):
    
        ## WHITE NOISES GENERATION
        
        wnoise   = the_white_noise(N)
        
        
        ## KRAVCHUK SPECTROGRAM ANALYSIS
        
        # compute the transform from Spin Spherical Harmonics
        Kn, thetas, phis   = the_ssht_transform(wnoise)
        
        # find the zeros
        znt, znp           = the_zeros(Kn,N,thetas,phis)

        
        ## SPATIAL STATISTICS ON ZEROS
        
        # explored range of r
        rs       = np.linspace(0,2*np.pi/np.sqrt(N),10**4)
        
        # compute the F-function
        Fn,_     = the_F_statistics(znt,znp,rs)
        
        # compute Ripley's K function
        Rn,_     = the_K_statistics(znt,znp)
        
        
        ## SAVE THE SPECTROGRAM AND THE SPATIAL FUNCTIONALS
        
        mdict    = {"N" : N, "wnoise" : wnoise, "Kn" : Kn, "znt" : znt, "znp" : znp, "rs" : rs, "Fn" : Fn, "Rn" : Rn}

        savemat(folder+'/noise_'+str(n)+'.mat',mdict)

        print('white noise '+str(n+1)+'/'+str(m)+' done.')
        
    return m, folder
