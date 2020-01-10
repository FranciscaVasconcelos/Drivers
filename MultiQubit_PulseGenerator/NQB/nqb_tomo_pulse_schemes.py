"""
 ----------------------------------------------------
 - Common pulse schemes for 1 and 2 qubit tomo      -
 ----------------------------------------------------

Functions and definitions used for determining which pulse scheme was used
to do the tomography and convert those pulses into corresponding Pauli matrices
"""

import numpy as np
import itertools

# Initialize Pauli matrices in numpy format:
I = np.matrix([[1, 0], [0, 1]])
sz = np.matrix([[1, 0], [0, -1]])
sx = np.matrix([[0, 1], [1, 0]])
sy = np.matrix([[0, -1j], [1j, 0]])

# Create dictionary used for mapping post-pulsing into axis
dictPulseToPauli = {'I': sz, 'X2p': sy, 'Y2m': sx, 'X2m': -sy, 'Y2p': -sx, 'i':I}



def NQubit_Pulse(n):
    """Function returning list of strings correponding to pulses used for
     two qubit 9 pulse tomo in Labber

    Returns
    -------
    PulseScheme_2QB_9Pulses (list):
        List of strings corresponding to pulses used in labber for the 9 pulse
        two qubit state tomography

    """

    PulseScheme_NQB = []
    for i in itertools.combinations(['Y2m-', 'X2p-', 'I-'], n):
        PulseScheme_NQB.append(i[:-1])
    return PulseScheme_NQB

    

def NQubit_Sim_Pulse(n):
    """Function returning list of strings correponding to pulses used for
     two qubit 9 pulse tomo in Labber

    Returns
    -------
    PulseScheme_2QB_9Pulses (list):
        List of strings corresponding to pulses used in labber for the 9 pulse
        two qubit state tomography

    """

    PulseScheme_NQB = []
    for i in list(itertools.product(['Y2m-', 'X2p-', 'I-', 'i-'], repeat=n)):
        i = ''.join(i)
        if i != "i-"*n: # exclude all identity case
            #print(i[:-1])
            PulseScheme_NQB.append(i[:-1])
    #print("Pulse Scheme: ",PulseScheme_NQB)
    return PulseScheme_NQB


# computes tensor product of given list of input
# order of tensors is order in list
def tensor_product(vals):
    #print(vals)
    tensor = vals[0]
    idx = 1

    #print(type(vals))

    while idx < len(vals):
        #print(vals[idx][0])
        #print(tensor)
        tensor = np.kron(np.matrix(tensor), np.matrix(vals[idx]))
        idx += 1

    return np.matrix(tensor)


def PulseToPauli_NQB(n, pulses):
    """Similar to PulseToPauli_1QB, except for 3-qubit pulsing for three-qubit
    state tomography. Returns 8x8 matrix corresponding to tensor product of
    pauli matrices

    Ex:
    Postpulse both qubits with -pi/2 on Y (Y2m-Y2m): corresponding 8x8 matrix
    is (sigma_x tensor sigma_x).

    Parameters
    ----------
    pulses : list of strings
        List containing names of pulses

    Note:
    Assumes standard Labber naming convention for pulses: 'Pulse1-Pulse2' where
    'PulseX' can be any of {I, X2p, X2m, Y2p, Y2m}

    Returns
    -------
    type
        Returns array of 4x4 matrices corresponding to outer products of Pauli
        matrices

    """

    tensor = 1

    Paulis_NQB = []

    for Pulse in pulses:
        comps = Pulse.split('-')
        matrices = []
        for c in comps:
            matrices.append(dictPulseToPauli[c])
        Paulis_NQB.append(tensor_product(matrices))
    #print('-------------------')

    #print(Paulis_NQB)

    return Paulis_NQB
