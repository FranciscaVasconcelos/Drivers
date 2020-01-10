import numpy as np
import matplotlib.pyplot as plt

def getTsFromCholesky(T):
	'''Function that creates list of t-parameters from cholesky T matrix.'''
	t = []
	for r in range(T.shape[0]):
		for c in range(T.shape[1]):
			if r == c:
				t.append(T[r][c].real)
			elif r>c:
				if T[r][c] == 0:
					t.append(0)
					t.append(0)
				elif isinstance(T[r][c], complex):
					t.append(T[r][c].real)
					t.append(T[r][c].imag)
				else:
					t.append(T[r][c].real)
					t.append(0)
				'''if c == 0:
					temp_dict[i] = [T[r][c].real]
				else:
					if isinstance(T[r][c], complex):
						temp_dict[i] = [T[r][c].real, T[r][c].imag]
					else: 
						temp_dict[i] = [T[r][c]]
				i+=1

	for key, val in sorted(temp_dict.items()):
		t = t + val'''
	print(len(t))
	return t

def getRhofromt_NQB(n,nlvl,t):
	"""Function that takes guess of Cholesky decomposition of density matrix,
	and returns a 8x8 density matrix

	The Cholesky decomposition for three-qubit density matrix is:
	rho = (T.dag() * T)/(Tr(T.dag() * T))
	where T is upper triangular matrix of the form
	T = [[  t_0       ,      0     ,    0      ,  0 ],
		 [  t_4+it_5  ,     t_1    ,    0      ,  0],
		 [  t_10+it_11,  t_6+it_7  ,   t_2     ,  0],
		 [  t_14+it_15, t_12+it_13 , t_8+it_9  , t_3] ...
		 ...]
	Parameters
	----------
	t : array
		length = 64, containing (real) values of t's

	Returns
	-------
	rho_t: array
		4x4 density matrix

	"""
	'''
	Function that takes a t (len=64) and generates the 8x8 density matrix from
	Cholesky decomp:rho = (T.dag() * T)/(Tr(T.dag() * T))

	'''
	'''T = np.matrix([[t[0],                 0,              0,              0],
				   [t[1] + 1j*t[2],     t[3],             0,              0],
				   [t[4] + 1j*t[5],  t[6] + 1j*t[7],    t[8],           0],
				   [t[9] + 1j*t[10], t[11] + 1j*t[12], t[13] + 1j*t[14], t[15]] ...])
	'''
	T = []
	k=0
	for r in range(nlvl**n):
		row = []
		for c in range(nlvl**n):
			if c > r:
				row.append(0)
			elif c == r:
				row.append(t[k])
				k+=1
			else:
				row.append(t[k]+1j*t[k+1])
				k+=2
		T.append(row)

	T=np.array(T)
	plt.matshow(T.real, cmap='hot', interpolation='nearest')
	plt.title("GENERATED CHOLESKY")
	plt.show()

	print(T.real)
	print("*****************************************")


	norm = np.dot(T, T.transpose().conj()).trace()
	print(norm)
	plt.matshow(T.conjugate().transpose().real, cmap='hot', interpolation='nearest')
	plt.title("CONJUGATE CHOLESKY")
	plt.show()
	rho_t = np.dot(T, T.transpose().conj())/norm
	return np.array(rho_t)

def idealFockRho(phase = 0):
    rho = np.zeros((4,4),dtype=complex)
    rho[1,1] = .5
    rho[2,2] = .5
    rho[2,1] = .5*np.exp(1j*phase)
    rho[1,2] = .5*np.exp(-1j*phase)
    return rho

def idealNOONRho(phase = 0):
	rho = np.zeros((9,9),dtype=complex)
	#print( 1e-1*np.eye(9))
	rho[2,2] = .5
	rho[6,6] = .5
	rho[2,6] = .5*np.exp(1j*phase)
	rho[6,2] = .5*np.exp(-1j*phase)
	return rho

#print(idealNOONRho().shape)
#print(np.eye(9).shape)

density_matrix = np.eye(9)
#density_matrix = idealNOONRho()
#density_matrix = idealFockRho()
plt.matshow(density_matrix.real, cmap='hot', interpolation='nearest')
plt.title("ORIG MATRIX")
plt.show()

try:
	t_guess = getTsFromCholesky(np.linalg.cholesky(density_matrix))
	plt.matshow(np.linalg.cholesky(density_matrix), cmap='hot', interpolation='nearest')
	plt.title("TRUE CHOLESKY")
	plt.show()
	print(np.linalg.cholesky(density_matrix))
except Exception as e:
	print(e)
	t_guess = getTsFromCholesky(np.linalg.cholesky(density_matrix+ 1e-14*np.eye(np.shape(density_matrix)[0])))
	plt.matshow(np.linalg.cholesky(density_matrix+ 1e-14*np.eye(np.shape(density_matrix)[0])).real, cmap='hot', interpolation='nearest')
	plt.title("TRUE CHOLESKY")
	plt.show()
	print(np.linalg.cholesky(density_matrix+ 1e-14*np.eye(np.shape(density_matrix)[0])).real)
	print("***************************************")


rho = np.array(getRhofromt_NQB(2,3,t_guess).real)
print("___________________________")
print(density_matrix)
print("___________________________")
print(rho)


plt.matshow(rho, cmap='hot', interpolation='nearest')
plt.title("GENERATED MATRIX")
plt.show()

print(rho)