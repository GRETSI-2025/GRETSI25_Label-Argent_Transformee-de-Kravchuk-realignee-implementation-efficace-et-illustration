# LOAD NECESSARY PYTHON LIBRARIES

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
import scipy.signal
import time


# IMPORT USEFUL FUNCTIONS

from scipy.io                        import wavfile


# LOAD THE FUNCTIONS OF THE KRAVCHUK TOOLBOX

from elementary_signals              import the_chirp, the_dirac, the_sine, the_white_noise, the_noisy_chirp, the_noisy_dirac, the_noisy_sine, the_coherent_state, the_noisy_coherent_state
from elementary_signals              import display_signal
from pkravchuk_transform             import the_ptransform
from ssht_transform                  import the_ssht_transform, the_new_transform, the_inverse_transform, rotate_signal
from kravchuk_display                import planar_display, spherical_display
from stft_transform                  import the_stft_transform, stft_display


# HANDLE APPEARANCE OF THE PLOTS

from mpl_toolkits.axes_grid1         import make_axes_locatable

mpl.rcParams['xtick.labelsize'] = 20;
mpl.rcParams['ytick.labelsize'] = 20;
mpl.rcParams['axes.titlesize'] = 20;
plt.rc('axes', labelsize=22.5);
plt.rc('legend', fontsize=20);
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
mpl.rcParams['font.family'] = 'roman'
