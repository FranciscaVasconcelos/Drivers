import sys
import os
import numpy as np
#import qutip as qt  # useful for calculating fidelities etc.
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt  # Useful for plotting
from matplotlib import gridspec as gridspec
import itertools
import random
import time

# Add directory to modules to path:
path_to_bfcmodules = '../bfcmodules/'
sys.path.append(path_to_bfcmodules)


import tomography_functions as tf
import pulse_schemes as ps 
import plotting_functions as pf

# THIS IS JUST SIMULATED DATA GENERATION

N = 3
	
s0 = np.array([[1],[0]])
s1 = np.array([[0],[1]])


sn0 = tf.tensor_product([s0]*N)
sn1 = tf.tensor_product([s1]*N)

state = (sn0+sn1)/np.sqrt(2)

state = np.outer(state, state)

density = np.zeros((2**N,2**N)).astype(complex)

pulse_scheme = ps.NQubit_Sim_Pulse(N)
paulis = ps.PulseToPauli_NQB(N,pulse_scheme)
probvectors = []

for p in paulis:
	expec = np.trace(state*p)
	mat = expec*p
	density += mat

density += np.identity(2**N)
density = density/(2**N)


for p in paulis:
	probvectors.append(np.trace(density.dot(p)).real)


# THIS IS THE CODE YOU WILL WANT TO USE

ts = tf.MLE_NQB_sim(N,np.array(probvectors), np.array(paulis), density, 
					tolerance = None, max_iter_num=20, random_initialization = False, 
					plot_convergence=True, verbose=True, directory_path="")

rho = tf.getRhofromt_NQB(N,ts.x) # The value of your density matrix!

# DENSITY MATRIX VISUALIZATION

fig = plt.figure()
gs = gridspec.GridSpec(1, 2)
ax1 = plt.subplot(gs[0, 0])
ax2 = plt.subplot(gs[0, 1])
#pf.plotRho3QB(rho, [ax1, ax2])

vals = itertools.product('01', repeat=N)
vals = [''.join(val) for val in vals]
x = [r'$|{0}\rangle$'.format(val) for val in vals]

pf.densityHeatMap(np.real(rho), x,x, "Real", ax1)
pf.densityHeatMap(np.imag(rho), x,x, "Imaginary", ax2)
plt.subplots_adjust(wspace=.4)
#plt.show()

fig.savefig('./'+str(N)+'_qubit.png')
plt.close()
