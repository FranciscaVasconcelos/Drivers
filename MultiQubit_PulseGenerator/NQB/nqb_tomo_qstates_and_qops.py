"""
 ----------------------------------------------------
 -  Useful shorthand definitions of states etc.     -
 ----------------------------------------------------

Just a list of useful single- and two-qubit matrices, rotation, states etc.

There's nothing conceptual in this file, just a lot of shorthands that
have proven useful
"""

import numpy as np
import qutip as qt

# initialize handy shortcuts using qutip
zero = qt.basis(2, 0)  # |0>
one = qt.basis(2, 1)  # |1>

# shorthand for Paulis
Z = qt.sigmaz()
X = qt.sigmax()
Y = qt.sigmay()
Id = qt.identity(2)

# shorthand for rotations
Y2p = qt.rotation(Y, np.pi / 2)
Y2m = qt.rotation(Y, -np.pi / 2)
Z2p = qt.rotation(Z, np.pi / 2)
Z2m = qt.rotation(Z, -np.pi / 2)
X2p = qt.rotation(X, np.pi / 2)
X2m = qt.rotation(X, -np.pi / 2)

# Initialize density matrices for fiducial singlequbit states
rho_0 = qt.ket2dm(zero)
rho_1 = qt.ket2dm(one)
rho_Xp = qt.ket2dm(Y2p * zero)
rho_Xm = qt.ket2dm(Y2m * zero)
rho_Yp = qt.ket2dm(X2m * zero)
rho_Ym = qt.ket2dm(X2p * zero)

# Initialize two qubit states
zerozero = qt.tensor(zero, zero)
zeroone = qt.tensor(zero, one)
onezero = qt.tensor(one, zero)
oneone = qt.tensor(one, one)

# Initialize density matrices for fiducial two-qubit states
# Return as numpy arrays
rho_00 = qt.ket2dm(zerozero).full()
rho_01 = qt.ket2dm(zeroone).full()
rho_10 = qt.ket2dm(onezero).full()
rho_11 = qt.ket2dm(oneone).full()

# Shorthand two-qubit operation you may be interested in:
U_2qb = qt.rotation(qt.tensor(Z, Z), np.pi / 2)

# Initialize ideal 2-qubit states
zerozero_U_2qb = U_2qb * qt.tensor(zero, zero)
zeroone_U_2qb = U_2qb * qt.tensor(zero, one)
onezero_U_2qb = U_2qb * qt.tensor(one, zero)
oneone_U_2qb = U_2qb * qt.tensor(one, one)

# Initialize density matrices for for fiducial two qubit states with your
# two-qubit of choice acting on it:
rho_00_U_2qb = qt.ket2dm(zerozero_U_2qb).full()
rho_01_U_2qb = qt.ket2dm(zeroone_U_2qb).full()
rho_10_U_2qb = qt.ket2dm(onezero_U_2qb).full()
rho_11_U_2qb = qt.ket2dm(oneone_U_2qb).full()
