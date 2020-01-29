'''
 -------------------------------------------------
 -      Suite of functions for tomography        -
 -------------------------------------------------

Generalized MLE code written by Francisca Vasconcelos (francisc@mit.edu),
based on 1-2QB code written by Morten Kjaergaard (mortenk@mit.edu)
with input from EQuS team and LL team.
'''
import numpy as np
from scipy.linalg import lu
import scipy
from scipy.optimize import minimize
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt 
from matplotlib.ticker import MaxNLocator
import MainDataHelpers.n_qubit_bharath.pulse_schemes as ps
import itertools
import time

from collections import Counter


converge_vals = {}
iter_density = []
iter_number = 0

TGLOBAL = 0
METHOD = 0
N = 0

def callbackIter(Xk):  
    global iter_density
    global TGLOBAL
    global METHOD
    global iter_number
    global N
    iter_density.append(Xk)
    cur_time = time.time()
    TGLOBAL = cur_time 
    iter_number += 1  
    if iter_number%10==0: print(iter_number)
       


def MLE_NQB_sim(n,  measurements,  paulis, variances, density_matrix, 
                tolerance, random_initialization, max_iter_num, nlvl=2,
                plot_convergence=False, verbose=False, directory_path=""):

    """Maximum likelihood estimation for n_qubit n_level qubit state tomography.

    For illustrative purposes, all examples are done with 3-qubit systems:

    Detailed description of the entire MLE process:
    The expectation value of 
        m_<III> 
        m_<IIA>
        m_<IAI>
        m_<IAA> 
        m_<AII>
        m_<AIA>
        m_<AAI>
        m_<AAA>
    for a 8x8 Matrix

    constructed as the tensor product of two paulis, can be found via:
    
    [m_<III>, m_<IIA>, m_<IAI>, m_<IAA>, m_<AII>, m_<AIA>, m_<AAI>, m_<AAA> ] 
        = betas^{-1}*[p_000^A, p_001^A, p_010^A, p_011^A, p_100^A, p_101^A, p_110^A, p_111^A]

    and the objective function L to be minimized to enfore PSD is given by:
    L = sum_{A}(m_<A> - Tr(A*rho_T) )**2

    where rho_T is the cholesky decomposition of a 8x8 matrix, and A is the
    tensor product of three pauli matrices. E.g. for the '9 pulse scheme'
    typically used in Labber:

    A =    [XXX = sx @ sx @ sx      YYY = sy @ sy @ sy      ZZZ = sz @ sz @ sz
            XXI = sx @ sx @ Id      YYI = sy @ sy @ Id      ZZI = sz @ sz @ Id
            XIX = sx @ Id @ sx      YIY = sy @ Id @ sy      ZIZ = sz @ Id @ sz
            XII = sx @ Id @ Id      YII = sy @ Id @ Id      ZII = sz @ Id @ Id
            IXX = Id @ sx @ sx      IYY = Id @ sy @ sy      IZZ = Id @ sz @ sz
            IXI = Id @ sx @ Id      IYI = Id @ sy @ Id      IZI = Id @ sz @ Id
            IIX = Id @ Id @ sx      IIY = Id @ Id @ sy      IIZ = Id @ Id @ sz
            
            XXY = sx @ sx @ sy      XXZ = sx @ sx @ sz      ZZY = sz @ sz @ sy
            XYX = sx @ sy @ sx      XZX = sx @ sz @ sx      ZYZ = sz @ sy @ sz
            XYY = sx @ sy @ sy      XZZ = sx @ sz @ sz      ZYY = sz @ sy @ sy
            YXX = sy @ sx @ sx      ZXX = sz @ sx @ sx      YZZ = sy @ sz @ sz
            YXY = sy @ sx @ sy      ZXZ = sz @ sx @ sz      YZY = sy @ sz @ sy
            YYX = sy @ sy @ sx      ZZX = sz @ sz @ sx      YYZ = sy @ sy @ sz

            XYI = sx @ sy @ Id      XZI = sx @ sz @ Id      ZYI = sz @ sy @ Id
            XIY = sx @ Id @ sy      XIZ = sx @ Id @ sz      ZIY = sz @ Id @ sy
            YIX = sy @ Id @ sx      ZIX = sz @ Id @ sx      YIZ = sy @ Id @ sz
            YXI = sy @ sx @ Id      ZXI = sz @ sx @ Id      YZI = sy @ sz @ Id
            IXY = Id @ sx @ sy      IXZ = Id @ sx @ sz      IZY = Id @ sz @ sy
            IYX = Id @ sy @ sx      IZX = Id @ sz @ sx      IYZ = Id @ sy @ sz

            XYZ = sx @ sy @ sz
            XZY = sx @ sz @ sy
            YXZ = sy @ sx @ sz
            YZX = sy @ sz @ sx
            ZXY = sz @ sx @ sy
            ZYX = sz @ sy @ sx]

    where @ corresponds to the tensor product. The corresponding pulse sequence
    to get the above list of A's is:
    Pulses = ['Y2m-Y2m-Y2m', ...]
    If you feed this function the pulse scheme in the format outlined above,
    it automatically generates the corresponding 4x4 matrices

    Parameters
    ----------
    n: int
        number of qubits

    probvectors : array
        Array of arrays of probabilities of the form 
            [p_000, p_001, p_010, p_011, p_100, p_101, p_110, p_111],
        for each three qubit state constructed as rho*(A @ B @ C), 
            where A, B, and C are three Pauli operators
    betas : array
        Array of arrays of betas for each qubit. Each beta corresponds to a
        matrix of the form
                    beta[0][0] = beta_I for Q1_|0>
                    beta[0][1] = beta_Z for Q1_|0>
                    beta[1][0] = beta_I for Q1_|1>
                    beta[1][1] = beta_Z for Q1_|1>
    pulse_scheme : Array of strings
        The pulse scheme used to generate the data. See the supplementary
        file 'pulse_schemes' for more descriptions
    verbose : bool
        If true will return not just minimal t values, but also more results
        from the minimizer


    PARAM RUNDOWN: 
            n = # qubits, 
            nlvl = # levels in system (default=2)
            measurements = list of expectation measurements,
            paulis = list of paulis (in same order as measurements),
            variances = 
            density_matrix = IDEAL form of output (what you are hoping for),
            tolerance = set to None if you don't care (can set specific value to increase speed),
            random_initialization = True (slower, use random initial guess), False (faster, use density matrix as initial guess)
            max_iter_num = 
            plot_convergence = True (make MLE convergence plots), False (don't make MLE convergence plots)
            verbose = verbose output (can leave as False),
            directory_path = if you want to save plots in a particular directory
    

    Returns
    -------
    Array of t's
        Array of 64 t's that minimize the cost function.

    """
    # just guessing all t parameters are 1/(nlvl^2)^n
    if random_initialization:
        t_guess = np.ones((nlvl**2)**n)/((nlvl)**2**n)  
    #set initial guess to be desired density matrix cholesky decomposition     
    else:
        try:
            t_guess = getTsFromCholesky(np.linalg.cholesky(density_matrix))
        except Exception as e: # if error, slightly perturb matrix
            t_guess = getTsFromCholesky(np.linalg.cholesky(density_matrix+ 1e-14*np.eye(np.shape(density_matrix)[0])))
    consts = ({'type': 'eq',
               'fun': lambda t: PSDconstraint_NQB(t)})

    # Now do some array juggling to make entry[0] in measuredEvals correspond
    # to the first 8x8 pauli matrix in Paulis. Since some pulses (e.g. 'I-I-I')
    # gives access to m_<ZZZ>, m_<ZZI>, m_<ZIZ>, m_<ZII>, m_<IZZ>, m_<IZI>, m_<IIZ> 
    # (which is then output from getEvalsfromProbs_3QB), 
    # the arrays need to be reshaped correctly:

    # Now input measured expectation values and paulis to the minimizer:
    if plot_convergence: # create convergence plots

        # paramters to keep track of during MLE process for plotting
        global iter_density # keep track of predicted params
        global TGLOBAL # keep track of wall clock time
        global iter_number # keep track of number of iterations
        global N # number of qubits

        N = n

        start_time = time.time()

        TGLOBAL = time.time()
        
        # normal MLE approach
        result = minimize(MLE_Functional_NQB,
                          t_guess,
                          args=(n,nlvl, measurements, paulis,variances),
                          constraints=consts,
                          options={'maxiter': max_iter_num, 'disp': False},
                          #tol=tolerance,
                          callback=callbackIter)

        iter_number = 0
        total_time = time.time()-start_time

        num_iterations = len(iter_density)

        print(n, "_", tolerance, "_", random_initialization, ": ", total_time)

        # convergence plot
        infidelities = []
        for i in iter_density:
            fidelity = np.trace(np.matmul(getRhofromt_NQB(n,nlvl,i) , density_matrix))
            infidelities.append(fidelity)
        f = open("./"+directory_path+"times.txt", "a")
        f.write("\n"+str(n)+" "+str(tolerance)+" "+str(random_initialization)+" "+str(num_iterations)+" "+str(total_time))
        fig = plt.figure()
        ax = fig.gca()
        plt.title("Tolerance="+str(tolerance)+" Random Init="+str(random_initialization), fontsize=10)
        plt.suptitle(str(n)+" Qubit Fidelity Convergence",fontsize=16, fontweight='bold')
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel("Iteration Number, K")
        plt.ylabel("Fidelity")
        plt.plot(range(1,len(infidelities)+1),infidelities)
        plt.plot(range(1,len(infidelities)+1),infidelities, 'ro')
        plt.xlim([0, max_iter_num])
        plt.ylim([0, 1.05])
        plt.plot([1 for n in range(max_iter_num)], ls='--', alpha=0.5, c='grey')
        fig.savefig("./"+directory_path+str(n)+"_"+str(tolerance)+"_"+str(random_initialization)+".png")
        plt.close()

        iter_density = []

    else:
        start_time = time.time()
        
        # normal MLE approach
        result = minimize(MLE_Functional_NQB,
                          t_guess,
                          args=(n,nlvl, measurements, paulis,variances),
                          constraints=consts,
                          options={'maxiter': max_iter_num, 'disp': False},
                          #tol=tolerance,
                          callback=callbackIter)


        total_time += time.time()-start_time

        print(n, "_", tolerance, "_", random_initialization, ": ", total_time)

    if verbose:  # in case you want more details than just the solution
        return result
    return result['x']


def MLE_Functional_NQB(t, n, measuredEvals, Paulis, variances, nlvl=2):
    """The functional to minimized for 2 qubit density matrix

    Parameters
    ----------
    t : array
        List of 16 floats, corresponding to each t_i
    measuredEvals : array
        Array of 15 floats, corresponding to measured expectation values of
        various Pauli operators
    Paulis : array
        Array of 15 4x4 Pauli matrices, corresponding to the measured
        expectation values

    Returns
    -------
    L : Float
        The number to minimized in the maximum likelihood minimizer (cost function)
    """
    r = getRhofromt_NQB(n,t, nlvl=nlvl)
    measuredEvals = np.array(measuredEvals)
    Paulis = np.array(Paulis)
    variances = np.array(variances)
    L = np.sum((np.abs(measuredEvals-np.array(
            [(r.dot(mat)).trace() for mat in Paulis]))/np.array(variances))**2)
    return L

def tensor_product(vals):
    """Function used for computing the tensor product of a given list of input tensors

    Parameters
    ----------
    vals : list
        list of tensors in order of tensor product

    Returns
    -------
    tensor
        tensor product output (np.array)

    """
    tensor = vals[0]
    idx = 1

    while idx < len(vals):
        tensor = np.kron(np.array(tensor), np.array(vals[idx]))
        idx += 1

    return tensor


def PSDconstraint_NQB(t):
    """Function used for enforcing positive semidefinite property of
     density matrix in cholesky decomposition

    Ensures that Tr(rho) = 1, which is equivalent to sum(t_i**2) = 1.

    Parameters
    ----------
    t : array
        t_i parameters that go into Cholesky decomposition

    Returns
    -------
    float
        Function that is 0 when PSD constraint satisfied

    """
    return np.array((t[:]**2).sum()-1)


def getRhofromt_NQB(n, nlvl, t):
    """Function that takes guess of Cholesky decomposition of density matrix,
    and returns a (nlvl^n)x(nlvl^n) density matrix

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
    n: int
        number of qubits
    nlvl: int
        number of levels in quantum system
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

    norm = np.dot(T, T.transpose().conj()).trace()
    rho_t = np.dot(T, T.transpose().conj())/norm
    return np.array(rho_t)

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

    return t

        






