import os
import sys
sys.path.append("../src")

import numpy                         as np
import matplotlib                    as mpl
import matplotlib.pyplot             as plt
import numpy.random                  as npr
import scipy.special                 as special
import statsmodels.stats.proportion  as stms
import pyssht                        as ssht
import h5py
import librosa
import librosa.display
import scipy.signal
import cmocean
import time

from mpl_toolkits.axes_grid1         import make_axes_locatable
from scipy.io                        import savemat, loadmat, wavfile

from chirp_signals          import the_chirp, the_gaussian, the_white_noise, the_noisy_chirp, the_noisy_gaussian, the_bird_call, the_noisy_bird
from chirp_signals          import the_sine, the_noisy_sine, display_signal, display_noisy_and_true, display_signal_and_zoom
from chirp_signals          import the_dirac, the_noisy_dirac, the_constant, the_noisy_constant
from kravchuk_transform     import the_transform, the_zeros
from kravchuk_display       import signal_display, planar_display, spherical_display
from pkravchuk_transform    import the_ptransform
from ssht_transform         import the_ssht_transform, the_new_transform, the_inverse_transform, rotate_signal
# from spherical_statistics   import the_distance, the_F_statistics, the_K_statistics, empirical_F, empirical_K
# from detection_test         import the_test_statistic, the_test, the_true_statistic, the_global_statistic
from stft_transform         import the_stft_transform, the_stft_zeros, stft_display
# from white_noises           import noise_samples
# from window_analysis        import analyze_signal, extract_results


mpl.rcParams['xtick.labelsize'] = 20;
mpl.rcParams['ytick.labelsize'] = 20;
mpl.rcParams['axes.titlesize'] = 20;
plt.rc('axes', labelsize=22.5);
plt.rc('legend', fontsize=20);
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
mpl.rcParams['font.family'] = 'roman'

sapin  = (0.0353, 0.3216, 0.1569)
sclair = (0.7588, 0.8304, 0.7892)
carmin = (0.7294,0.0392,0.0392)
bleu   = (0.2, 0.2, 0.7020)
bclair = (0.6235, 0.6235, 0.8980)
rose   = (0.56, 0.004, 0.32)
marron = (0.51, 0.3, 0.11)
jaune  = (0.93, 0.69, 0.13)
tendre = (0.106, 0.463, 0.827)
tclair = (0.776, 0.866, 0.957)


my_cmap = cmocean.cm.cmap_d['deep']
c       = np.zeros((5,3))
rgba    = my_cmap(0.999)
c[0,:]  = rgba[:3]
rgba    = my_cmap(0.8)
c[1,:]  = rgba[:3]
rgba    = my_cmap(0.6)
c[2,:]  = rgba[:3]
rgba    = my_cmap(0.4)
c[3,:]  = rgba[:3]
rgba    = my_cmap(0.2)
c[4,:]  = rgba[:3]


d         = np.zeros((7,3))
d[0,:]    = np.array([0,0,0])
d[1:-1,:] = c
rgba      = my_cmap(0.1)
d[-1,:]   = rgba[:3]