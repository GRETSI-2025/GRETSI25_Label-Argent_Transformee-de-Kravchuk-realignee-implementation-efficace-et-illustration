import numpy              as np
import matplotlib.pyplot  as plt
import os
import cmocean

from scipy.io             import savemat, loadmat
from spherical_statistics import empirical_F, empirical_K
from chirp_signals        import display_signal

my_cmap = cmocean.cm.cmap_d['deep']
rgba    = my_cmap(0.6)
my_blue = rgba[:3]

def the_test_statistic(S,S0,the_norm=2,r2 = int(1e4)):
    
    """
    test statistics with possible choice of the norm
    """
    
    # range of rs of interest
    ind_rs = np.arange(0,r2)
    
    # test statistics
    T = np.linalg.norm(S[ind_rs]-S0[ind_rs],ord=the_norm)
    
    return T
    
    
def the_test(signal,alpha=0.05,m=199,folder='samples',functional='F'):
    
    """
    summary statistics used to perform the test
    """
    
    #signal: observed data to be tested
    #alpha: confidence level of the test
    #folder: folder in which the white 
    #functional: either 'F' for the empty space function or 'K' for Ripley's K function
    
    k = int(alpha*(m+1))
    
    # Kravchuk analysis and functional statistics of the data
    if functional     == 'F':
        S,_           = empirical_F(signal)
    else:
        S,_           = empirical_K(signal)
    nr        = len(S)
    S_m       = np.zeros((m+1,nr))
    S_m[-1,:] = S
    
    # load the functional statistics of the noise samples
    for n in np.arange(m):
        ndict         = loadmat('../'+folder+'/noise_'+str(n)+'.mat',squeeze_me=True);
        if functional == 'F':
            S_m[n,:]  = ndict["Fn"]
        else:
            S_m[n,:]  = ndict["Rn"]
    
    # compute the summary statistics
    S0                = np.sum(S_m, axis = 0)/(m+1)
    T_m               = np.zeros(m)
    for i in np.arange(m):
            T_m[i]    = the_test_statistic(S_m[i,:],S0)
    T_m.sort()
    t_exp             = the_test_statistic(S,S0)
    
    if t_exp >= T_m[-k]:
        print('The null hypothesis is rejected with confidence '+str(1 - alpha)+'. A signal is detected.')
    else:
        print('The null hypothesis cannot be rejected with confidence '+str(1 - alpha)+'. No signal detected.')
    
    print(' ')
    print('--------------------------------------------------------------')
    print(' ')
    print('Threshold of rejection: %.2f' % (T_m[-k]))
    print('Value of the summary statistics: %.2f' % t_exp)
    
    return t_exp


def the_true_statistic(folder, Nreal = 10, functional = 'F', disp = False):

    """
    compute the ground truth for the chosen spatial statistic over available realizations of white noise
    """

    #folder: where to find the analyzed white Gaussian noises
    #Nreal: number of batches of realizations of white noise to consider
    #functional: either F for empty space function or R for Ripley's K function

    # Load one first example
    wdict           = loadmat(folder + '/noise_0.mat',squeeze_me=True)
    
    # Intialize the ground truth chosen spatial functional of the zeros of the spherical GAF
    Ntot            = 0
    if functional == 'F':
        Lts         = np.zeros((200 * Nreal, len(wdict["Fn"])))
        Lts[Ntot,:] = wdict["Fn"]
    else:
         Lts         = np.zeros((200*Nreal,len(wdict["Rn"])))
         Lts[Ntot,:] = wdict["Rn"]
    Ntot            += 1
        
    # All realizations available in the first batch
    for n in range(1,200):
        wdict           = loadmat(folder + '/noise_' + str(n) + '.mat',squeeze_me=True)
        if functional   == 'F':
            Lts[Ntot,:] = wdict["Fn"]
        else:
            Lts[Ntot,:] = wdict["Rn"]
        Ntot           += 1
    print('Batch 1/' + str(Nreal) + ' done.')
    
    # Additional realizations
    for real in range(1,Nreal):
        for n in range(200):
            wdict           = loadmat(folder + '/real-' + str(real) + '/noise_' + str(n) + '.mat',squeeze_me=True)
            if functional   == 'F':
                Lts[Ntot,:] = wdict["Fn"]
            else:
                Lts[Ntot,:] = wdict["Rn"]
            Ntot           += 1
        print('Batch ' + str(real + 1) + '/' + str(Nreal) + ' done.')
        
    # Compute mean and confidence region on all considered samples
    Lt              = np.mean(Lts,axis = 0)
    sLt             = 1.96 * np.std(Lts,axis = 0) / np.sqrt(Ntot)

    print('Total number of samples considered: ' + str(Ntot) + '.')

    # Display the empirical estimate of the chosen spatial functonal
    if disp:
        plt.subplots(figsize=(6, 3))
        if functional == 'F':
            plt.plot(wdict["rs"],Lt,linewidth = 2, color = my_blue, label  = r'$\overline{F}$')
        else:
            plt.plot(wdict["rs"],Lt,linewidth = 2, color = my_blue, label  = r'$\overline{K}$')
        plt.grid()
        plt.xlabel('$r$: geodesic distance')
        if functional == 'F':
            plt.title('Empirical estimate of empty space function')
        else:
            plt.title('Empirical estimate of Ripley''s $K$ function')
        plt.legend(loc = 4);

    return Lt, sLt

def the_global_statistic(snr, real, duration = 0.05, Lt = [], functional = 'F',test = 'sum', disp = False, root = '/Volumes/LaCie/ownCloud/Python/ssht-kravchuk-transform'):
    
    """
    global summary statistics obtained by scanning the signal
    """

    #snr: signal-to-noise ratio of the signal on which to perform the test
    #real: realization of the signal
    #duration: proportion of the observation window spanned by the chirp signal of interest
    #Lt: ground truth for the chosen spatial functional
    #functional: either F for empty space function or R for Ripley's K function
    #test: sum or max for the two different summary statistics
    #disp: if true display the signal
    #root: location of the simulation results

    # if not provided load the ground truth for the chosen spatial functional
    if len(Lt) == 0:
        nfolder    = root + '/data/noises'
        Nreal      = 20
        file       = nfolder + '/the_true_statistic_' + functional + '_' + str(Nreal) + '.mat'
        if os.path.isfile(file):
            Ldict      = loadmat(file, squeeze_me = True)
            Lt         = Ldict["Lt"]
        else:
            Lt         = the_true_statistic(nfolder, Nreal, functional)
            Ldict      = {"Lt" : Lt, "functional" : functional, "Nreal" : Nreal}
            savemat(folder + '/the_true_statistic_' + functional + '_' + str(Nreal) + '.mat',Ldict)
        
    # load the signal
    snrs   = str(int(100*snr))
    if duration == 0.05:
        sdict  = loadmat(root + '/data/signals/snr_'+snrs.rjust(3,'0')+'/synthetic_bird_call_'+str(real)+'.mat',squeeze_me=True)
    else:
        durs   = str(int(100*duration))
        sdict  = loadmat(root + '/data/signals/duration_'+durs.rjust(3,'0')+'/snr_'+snrs.rjust(3,'0')+'/synthetic_bird_call_'+str(real)+'.mat',squeeze_me=True)
    time_t = sdict["time_t"]
    signal = sdict["signal"]

    # indicate where to find the sliding window analysis results
    if duration == 0.05:
        folder         = root + '/data/signals/snr_'+snrs.rjust(3,'0')+'/real_'+str(real)
    else:
        folder         = root + '/data/signals/duration_'+durs.rjust(3,'0')+'/snr_'+snrs.rjust(3,'0')+'/real_'+str(real)
    name           = 'wsynthetic_bird_call_'+str(real)

    # display it
    if disp:
        time_t = display_signal(signal,time_t)

    # characteristics of the overlapping windows
    win        = 1025 
    s          = 32

    # initialize indices and lists
    Nwin    = 0
    start   = 0
    S_exp   = []
    
    while start + win < len(signal):
    
        # load the analysis of the Nwin window of the signal
        sdict      = loadmat(folder+'/'+name+'_'+str(Nwin)+'.mat',squeeze_me = True)
        if functional == 'F':
            Lb     = sdict["Fb"]
        else:
            Lb     = sdict["Rb"] 
    
        # compute the local L2 difference between empirical functional and ground truth
        s_exp      = np.linalg.norm(Lb-Lt)

        # store it
        S_exp.append(s_exp)
        
        # move to next window
        start    += s
        Nwin     += 1

    if test == 'sum':
        S_tot = np.sum( np.array(S_exp) )
    elif test == 'max':
        S_tot = np.max( np.array(S_exp))
    else:
        error('Summary stastistic not implemented.')

    return S_tot, S_exp
        
        