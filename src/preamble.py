# LOAD NECESSARY PYTHON LIBRARIES

import sys

sys.path.append("../src")

import matplotlib as mpl
import matplotlib.pyplot as plt


# IMPORT USEFUL FUNCTIONS


# LOAD THE FUNCTIONS OF THE KRAVCHUK TOOLBOX


# HANDLE APPEARANCE OF THE PLOTS


mpl.rcParams["xtick.labelsize"] = 20
mpl.rcParams["ytick.labelsize"] = 20
mpl.rcParams["axes.titlesize"] = 20
plt.rc("axes", labelsize=22.5)
plt.rc("legend", fontsize=20)
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
mpl.rcParams["font.family"] = "roman"
