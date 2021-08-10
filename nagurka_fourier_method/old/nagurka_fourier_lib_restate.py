from icecream import ic



from math import pi
import numpy as np

def P_restate_coeffs(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
    '''
    as opposed to the other formulation, here the coefficients make real sense
    '''
    m = np.arange(1, M+1)
    p = np.zeros((6))

    p0 = X_t0 - np.sum(a)#hey, why doesn't b factor in here?
    pf = X_tf - np.sum(a)# because sin at the edges is 0, so that term drops out

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


    '''
    Wait, there's even another restatement using alpha and beta coefficients in the fourier summation

    So, 4 papers and ~three conference proceedings, each using a different statement of the same technique.

    1988 American Control Conference - fifth order polynomial,
    alphas and beta coefficients on the sine terms

    1988 ASME publication presented at Winter Annual Meeting -
    third order polynomial, four boundary conditions,
    alpha and beta coefficients on the sine terms

    there's also a "multi segment" verson
    '''


# tests for this restatement file
# from nagurka_fourier_lib_restate import P_restate, P_restate_coeffs
# def test_P_restate_coeffs():
#     t_f = 2.0
#
#     a = np.array([-1.95643e-3])
#     b = np.array([1.442172e-3])
#     M = 1
#
#     t = np.linspace(epsilon, t_f, 100)
#
#     virus = default_virus(t)
#     host_cell = default_host_cell(t)
#
#     X_t0 = 1.0 # is it possible that this is overconstrained?
#     X_tf = 2.0
#     d_X_t0 = 0
#     d_X_tf = 3.0
#     d_d_X_t0 = 0
#     d_d_X_tf = -3.0
#
#     p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = P_restate_coeffs(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#
#     check.almost_equal(p0, X_t0 - L_(0.0,a,b,M,t_f))
#     check.almost_equal(pf, X_tf - L_(t_f,a,b,M,t_f))
#
#     check.almost_equal(d_p0, d_X_t0 - d_L_(epsilon,a,b,M,t_f), 1e-4)
#     check.almost_equal(d_pf, d_X_tf - d_L_(t_f,a,b,M,t_f), 1e-4)
#     check.almost_equal(d_d_p0, d_X_t0 - d_d_L_(epsilon,a,b,M,t_f), 1e-4)
#     check.almost_equal(d_d_pf, d_X_tf - d_d_L_(t_f,a,b,M,t_f), 1e-4)


# def test_compare_P():
#
#     a = np.array([1.0])
#     b = np.array([-0.3])
#     M = 1
#
#     t = np.linspace(epsilon, t_f, 100)
#
#     virus = default_virus(t)
#     host_cell = default_host_cell(t)
#
#     X_t0 = 0.0
#     X_tf = 5.0
#     d_X_t0 = 0.0
#     d_X_tf = 2.0
#     d_d_X_t0 = 6.1025137
#     d_d_X_tf = -3.4798053
#
#     P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#
