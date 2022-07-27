from ase.io import read
import numpy as np
import os,sys
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import pickle as pkl
from matplotlib.colors import Normalize
from matplotlib import cm

home='/Users/geokast/Papers/PdH_CO2/'
sys.path.append(home+'scripts') #Adapt this path your own
from intermediates_dict import *
from general_tools import quad_fun,lin_fun,get_reference_energies


plt.rcParams["figure.figsize"] = (7,5)
plt.rcParams["font.family"] = "Times New Roman"
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (7,5)
markersize=4

pj=os.path.join

pkldir=home+'results/parsed_all_data.pkl'
pkldir2=home+'results/parsed_most_stable_data.pkl'
rundir = pj(home,'data')
