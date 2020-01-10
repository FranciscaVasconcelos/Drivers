'''
 ----------------------------------------------------
 - Plotting functions for 1- and 2-qubit tomograpy  -
 ----------------------------------------------------

Mainly written by Morten Kjaergaard (mortenk@mit.edu)
with a ton of input from EQuS team and LL team.
'''
import numpy as np
import matplotlib.pyplot as plt

def densityHeatMap(array, x_label, y_label, title, ax):

    im = ax.imshow(array)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(x_label)))
    ax.set_yticks(np.arange(len(y_label)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(x_label)
    ax.set_yticklabels(y_label)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(y_label)):
        for j in range(len(x_label)):
            #print(array[i, j])
            text = ax.text(j, i, round(array[i, j], 2), ha="center", va="center", color="w",fontsize=6)


    ax.set_title(title)
    


def addOriginPointer(ax):
    """Adds marker at origin of axis instace

    Parameters
    ----------
    ax : axis instance
        axis on which to add origin maker

    """
    ax.axhline(0, c='#606060', linestyle='dashed', linewidth=0.2)
    ax.axvline(0, c='#606060', linestyle='dashed', linewidth=0.2)
    ax.scatter(0, 0, 25, c='#606060')


def generateTomoAxis_1QB(ax):
    """Modifies 3d axis instance for nice single qubit density matrix plotting

    Changes x- and y-tick locations, labels and limits

    Parameters
    ----------
    ax : 3d axis instance
        3D axis instance, e.g. ax = fig.add_subplot(12, projection='3d')
    """
    ticks = np.arange(0.5, 2, 1)
    stateticks = [r'$|0\rangle$', r'$|1\rangle$']
    ax.set_xticks(ticks)
    ax.set_xticklabels(stateticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(stateticks)
    ax.set_xlim([0.15, 1.85])
    ax.set_ylim([0.15, 1.85])
    ax.set_zlim([-1.1, 1.1])


def generateTomoAxis_2QB(ax):
    """Modifies 3d axis instance for nice two qubit density matrix plotting

    Parameters
    ----------
    ax : 3d axis instance
        3D axis instance, e.g. ax = fig.add_subplot(12, projection='3d')
    """
    ticks = np.arange(0.5, 4, 1)
    stateticks = [r'$|00\rangle$', r'$|01\rangle$',
                  r'$|10\rangle$', r'$|11\rangle$']
    ax.set_xticks(ticks)
    ax.set_xticklabels(stateticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(stateticks)
    ax.set_xlim([0.15, 3.85])
    ax.set_ylim([0.15, 3.85])
    ax.set_zlim([-1.1, 1.1])

def generateTomoAxis_3QB(ax):
    """Modifies 3d axis instance for nice three qubit density matrix plotting

    Parameters
    ----------
    ax : 3d axis instance
        3D axis instance, e.g. ax = fig.add_subplot(12, projection='3d')
    """
    ticks = np.arange(0.5, 8, 1)
    stateticks = [r'$|000\rangle$', r'$|001\rangle$', r'$|010\rangle$', r'$|011\rangle$',
                  r'$|100\rangle$', r'$|101\rangle$', r'$|110\rangle$', r'$|111\rangle$']
    ax.set_xticks(ticks)
    ax.set_xticklabels(stateticks)
    ax.set_yticks(ticks)
    ax.set_yticklabels(stateticks)
    ax.set_xlim([0.15, 7.85])
    ax.set_ylim([0.15, 7.85])
    ax.set_zlim([-1.1, 1.1])


def generateSphereAxis(sphere):
    """Overwrite standard axis of Bloch() to something neater

    Changes label convention, the sphere and frame alpha and gets rid of
    dangling Z label

    Parameters
    ----------
    sphere : Bloch() instance from QuTiP
        sphere for which to update label.

    """
    sphere.set_label_convention("xyz")
    sphere.sphere_alpha = 0
    sphere.frame_alpha = 0.2
    sphere.zlabel = ['', '']


def addVecsToSphere(sphere, vector):
    """ Adds vector to qutip bloch sphere instance

    Adds vector in black, and the x-, y- and z-axis projections in C0,C1 and C2
    colors

    Parameters
    ----------
    sphere : Bloch() instance from QuTiP
        Sphere on which to add vector`.
    vector : array
        Array with structure [<X>,<Y>,<Z>].

    """
    norm = np.linalg.norm(vector)  # calculate length of vector
    sphere.add_vectors(vector)  # Add initial polarization
    # Add explicit x,y,z vector
    sphere.add_vectors([vector[0], 0, 0])
    sphere.add_vectors([0, vector[1], 0])
    sphere.add_vectors([0, 0, vector[2]])
    # Print |p| as z-axis  name:
    sphere.zlabel = [r'$|\vec p|$ ={:2.2f}'.format(norm), '']
    sphere.vector_color = ['k', 'C0', 'C1', 'C2']


def plot_SS(ax, Data_0, Data_1, colors=None, QBnum=None):
    """Plots single shot data in axis instance ax

    Parameters
    ----------
    ax : axis instance
        Axis on which to plot the single shot data
    Data_0 : array
        Single shot data for |0>. Structure:
            [[Re(x0),Re(x1),...],[Im(x0),Im(x1),...]]
    Data_1 : array
        Single shot data for |1>
    """
    addOriginPointer(ax)
    if colors is not None:
        ax.scatter(Data_0[:, 0], Data_0[:, 1], 2, c=colors[0],
                   label=r"$|0\rangle$")
        ax.scatter(Data_1[:, 0], Data_1[:, 1], 2, c=colors[1],
                   label=r"$|1\rangle$")
    else:
        ax.scatter(Data_0[:, 0], Data_0[:, 1], 2, c="C0",
                   label=r"$|0\rangle$")
        ax.scatter(Data_1[:, 0], Data_1[:, 1], 2, c="C3",
                   label=r"$|1\rangle$")
    ax.set_xlabel("$I$ (mV)")
    ax.set_ylabel("$Q$ (mV)")
    ax.legend(loc=2)
    if QBnum is not None:
        ax.set_title('Single shot data for QB'+str(QBnum))
    else:
        ax.set_title('Single shot data')
    ax.tick_params(tickdir='in')


def plot_SVM(ax, SVM, Data_0, Data_1, colors=None, QBnum=None):
    """Plots a SVM used to discern optimal separating hyperplane between
    Data_0 and Data_1

    Parameters
    ----------
    ax : axis instance
        Axis instance on which to plot
    SVM : svm.SVC instance (from the sklearn package)
        trained SVC instance to be plotted
    Data_0 : array
        Data corresponding to |0>-state
    Data_1 : type
        Data corresponding to |1>-state

    """
    addOriginPointer(ax)
    if colors is not None:
        ax.scatter(Data_0[:, 0], Data_0[:, 1], 3, c=colors[0],
                   label=r"$|0\rangle$, training data")
        ax.scatter(Data_1[:, 0], Data_1[:, 1], 3, c=colors[1],
                   label=r"$|1\rangle$, training data")
    else:
        ax.scatter(Data_0[:, 0], Data_0[:, 1], 3, c="C0",
                   label=r"$|0\rangle$, training data")
        ax.scatter(Data_1[:, 0], Data_1[:, 1], 3, c="C3",
                   label=r"$|1\rangle$, training data")

    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    # create grid to evaluate model stored in SVC
    xx = np.linspace(xmin, xmax, 100)
    yy = np.linspace(ymin, ymax, 100)
    YY, XX = np.meshgrid(yy, xx)
    xy = np.vstack([XX.ravel(), YY.ravel()]).T  # annoying footwork
    Z = SVM.decision_function(xy).reshape(XX.shape)  # evaluate SVM on grid

    # plot decision boundary and margins
    ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5,
               linestyles=['--', '-', '--'], linewidths=[1]*3)
    # plot support vectors:
    ax.scatter(SVM.support_vectors_[:, 0], SVM.support_vectors_[:, 1],
               marker='o', s=1, color='white')
    # following line is just to make a nice legend. Very hacky
    ax.plot([], [], '.', marker='o', color='white', fillstyle='full',
            markeredgecolor='black', markeredgewidth=1.0, label='support')

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_xlabel("$I$ (mV)")
    ax.set_ylabel("$Q$ (mV)")
    ax.legend(loc=2)
    ax.set_title(r'Training SVM to discern 0 and 1')
    if QBnum is not None:
        ax.set_title(r'Training SVM to discern 0 and 1, for QB'+str(QBnum))


def plotBetaCalibration(ax, result_p0, result_p1, colors=None, QBnum=None):
    """Plots result of fitting cosines to beta calibration data

    Parameters
    ----------
    ax : axis instance
        Axis instance to plot results in
    result_p0 : lmfit.model
        Contains result of fitting to p0
    result_p1 : lmfit.model
        Contains result of fitting to p1

    """
    xData = result_p0.userkws['x']  # this is why lmfit is so damn awesome
    # Plot both p0 and p1 and write the result of fitting in the legend:
    ax.plot(xData, result_p0.best_fit, '-', color='k',
            label='fit, p0, [{:4.4f},{:4.4f}]'.format(
                result_p0.best_values['constant'],
                result_p0.best_values['amplitude']))
    ax.plot(xData, result_p1.best_fit, '-', color='gray',
            label='fit, p1, [{:4.4f},{:4.4f}]'.format(
                result_p1.best_values['constant'],
                result_p1.best_values['amplitude']))
    if colors is not None:
        ax.plot(xData, result_p0.data, 'o', label=r'p0', c=colors[0])
        ax.plot(xData, result_p1.data, 'o', label=r'p1', c=colors[1])
    else:
        ax.plot(xData, result_p0.data, 'o', label=r'p0', c='C0')
        ax.plot(xData, result_p1.data, 'o', label=r'p1', c='C3')

    if QBnum is not None:
        ax.set_title('Beta calibration for QB'+str(QBnum))
    else:
        ax.set_title('Beta calibration')

    ax.set_xlabel('Variable for Rabi pulse')
    ax.set_ylabel('Probabilites')
    ax.set_ylim([-0.01, 1.01])
    ax.legend(loc=1)

def plot_rho_1QB(rho, ax1, ax2):
    """'Skyscraper'-style plotting routine for 2x2 density matrix

    Plots real and imaginary part in two side-by-side panels
    Parameters
    ----------
    rho : numpy array
        2x2 numpy array to plot
    ax1 : 3d axis instance
        axis instance to plot real part of rho
    ax2 : type
        axis instance to plot real part of rho

    """
    x_pts, y_pts = np.meshgrid(np.arange(rho.shape[1]),
                               np.arange(rho.shape[0]))
    x_pts = x_pts.flatten()
    y_pts = y_pts.flatten()
    z_pts = rho.flatten()

    generateTomoAxis_1QB(ax1)
    generateTomoAxis_1QB(ax2)

    Re_rho_data = np.real(rho).flatten()
    Im_rho_data = np.imag(rho).flatten()
    # For aesthetics we offset points a tiny point to leave whitespace
    ax1.bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Re_rho_data, color='C0')
    ax2.bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Im_rho_data, color='red')

    ax1.set_title(r'Re($\rho$)')
    ax2.set_title(r'Im($\rho$)')
    plt.tight_layout()


def plotRho2QB(rho, two_ax):
    """'Skyscraper'-style plotting routine for 4x4 density matrix

    Plots real and imaginary part in two side-by-side panels
    Parameters
    ----------
    rho : numpy array
        4x4 numpy array to plot
    two_ax : Array of 2 '3d' axis instances
        First axis will be used for real part, 2nd for im. part
    """
    x_pts, y_pts = np.meshgrid(np.arange(rho.shape[1]),
                               np.arange(rho.shape[0]))
    x_pts = x_pts.flatten()
    y_pts = y_pts.flatten()
    z_pts = rho.flatten()

    generateTomoAxis_2QB(two_ax[0])
    generateTomoAxis_2QB(two_ax[1])

    Re_rho_data = np.real(rho).flatten()
    Im_rho_data = np.imag(rho).flatten()
    two_ax[0].bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Re_rho_data, color='C0')
    two_ax[1].bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Im_rho_data, color='red')
    two_ax[0].set_title(r'Re($\rho$)')
    two_ax[1].set_title(r'Im($\rho$)')

def plotRho3QB(rho, three_ax):
    """'Skyscraper'-style plotting routine for 8x8 density matrix

    Plots real and imaginary part in two side-by-side panels
    Parameters
    ----------
    rho : numpy array
        8x8 numpy array to plot
    three_ax : Array of 2 '3d' axis instances
        First axis will be used for real part, 2nd for im. part
    """
    x_pts, y_pts = np.meshgrid(np.arange(rho.shape[1]),
                               np.arange(rho.shape[0]))
    x_pts = x_pts.flatten()
    y_pts = y_pts.flatten()
    z_pts = rho.flatten()

    generateTomoAxis_3QB(three_ax[0])
    generateTomoAxis_3QB(three_ax[1])

    Re_rho_data = np.real(rho).flatten()
    Im_rho_data = np.imag(rho).flatten()
    three_ax[0].bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Re_rho_data, color='C0')
    three_ax[1].bar3d(x_pts + 0.1, y_pts + 0.1, np.zeros(len(z_pts)), 0.7, 0.7,
              Im_rho_data, color='red')
    three_ax[0].set_title(r'Re($\rho$)')
    three_ax[1].set_title(r'Im($\rho$)')
