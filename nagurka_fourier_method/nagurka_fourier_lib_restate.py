from icecream import ic

# def P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
#     '''
#     Eqs 22 to 27, Nagurka&Yen 1990
#     1988 has an alternate setup
#     '''
#
#
#     return p
#
#



from math import pi
import numpy as np


def P_restate_coeffs(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
    '''
    as opposed to the other formulation, here the coefficients make real sense
    '''
    m = np.arange(1, M+1)
    p = np.zeros((6))

    p0 = X_t0 - np.sum(a)
    pf = X_tf - np.sum(a)

    d_p0 = (d_X_t0 + 2.0 * np.sum(pi*m*b / t_f))
    d_pf = (d_X_tf + 2.0 * np.sum(pi*m*b / t_f))

    d_d_p0 = d_d_X_t0 + 4.0 * np.sum((pi**2.0)*(m**2.0)*a / (t_f**2.0))
    d_d_pf = d_d_X_tf + 4.0 * np.sum((pi**2.0)*(m**2.0)*a / (t_f**2.0))

    return p0, pf, d_p0, d_pf, d_d_p0, d_d_pf

def P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
    p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = P_restate_coeffs(t, X_t0, d_X_t0, d_d_X_t0, X_tf,
                                d_X_tf, d_d_X_tf, t_f, a, b, M)

    # how can I stop making these stupid mistakes?
    '''
    Eq 6 Nagurka and Yen and Benaroya 1987

    A Fourier-Based Method for the Suboptimal Control of Nonlinear Dynamical Systems
    '''
    P = t*0.0 # supports both float and array type
    tau = t/t_f

    P = ((-6.0*(p0-pf))-(3.0*(d_p0+d_pf)*(t_f))-(0.5*(d_d_p0 - d_d_pf)*(t_f**2.0)))*(tau**5.0)

    P += ((15.0*(p0-pf))+((8.0*d_p0+7.0*d_pf)*(t_f))+(0.5*((3/2.0)*d_d_p0 - d_d_pf)*(t_f**2.0)))*(tau**4.0)
    P += ((-10.0*(p0-pf))-((6.0*d_p0+4.0*d_pf)*(t_f))-(0.5*(3.0*d_d_p0 - d_d_pf)*(t_f**2.0)))*(tau**3)
    P += ((0.5*d_d_p0*(t_f**2.0)))*(tau**2)
    P += (d_p0*t_f)*(tau)
    P += p0
    return P
