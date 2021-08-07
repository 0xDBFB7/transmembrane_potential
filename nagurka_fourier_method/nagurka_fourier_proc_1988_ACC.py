from icecream import ic

'''

"A fourier based optimal control approach for structural systems",
Proceedings of the 1988 American Control Conference

eqs. 5 to 10.

'''



from math import pi
import numpy as np


def P_coefficients(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):

    m = np.arange(1, M+1)
    p = np.zeros((6))

    p0 = X_t0 - np.sum(a)
    pf = X_tf - np.sum(a)

    d_p0 = (d_X_t0 + 2.0 * np.sum(pi*m*b / t_f))
    d_pf = (d_X_tf + 2.0 * np.sum(pi*m*b / t_f))

    d_d_p0 = d_d_X_t0 + 4.0 * np.sum((pi**2.0)*(m**2.0)*a / (t_f**2.0))
    d_d_pf = d_d_X_tf + 4.0 * np.sum((pi**2.0)*(m**2.0)*a / (t_f**2.0))

    return p0, pf, d_p0, d_pf, d_d_p0, d_d_pf

def P_(t, p, M):
    p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = P_coefficients(t, X_t0, d_X_t0, d_d_X_t0, X_tf,
                                d_X_tf, d_d_X_tf, t_f, a, b, M)

    # how can I stop making these stupid mistakes?

    P = t*0.0 # supports both float and array type
    tau = t/t_f

    P =

    return P


def L_():
