"""
Test cases specifically applying the Yen parameterization to the transmembrane lib.
"""
'''

N.B.

Contains the actual virus optimization routine.

'''

from nagurka_fourier_lib import *

import sys
sys.path.app('../')
from transmembrane_lib import *


def X_to_


def test_transmembrane_compare():
    t_f = 0.001

    a = np.array([1.0, 4.0])
    b = np.array([-0.3,-2])
    M = 2

    t = np.linspace(epsilon, t_f, 100)

    virus = default_virus(t)
    default_host_cell = default_host_cell(t)


    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 4.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)



M = 1
