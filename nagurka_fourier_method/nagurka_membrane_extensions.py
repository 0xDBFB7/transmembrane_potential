
from nagurka_fourier_lib import *

import sys
sys.path.append('../')
from transmembrane_lib import *

from scipy.integrate import odeint
from scipy import signal

from scipy.optimize import minimize,basinhopping


def U_to_X(U, t, cell):
    yout = convolve_output(U,  cell, t[1]-t[0])
    # system = signal.lti(np.float64([cell.a_3,cell.a_2,cell.a_1]),np.float64([cell.b_3,cell.b_2,cell.b_1]))
    #
    # _, yout, _ = signal.lsim2(system, U=U, T=t, rtol=1e-11,atol=1e-11)#lsim2 ,
    return yout
