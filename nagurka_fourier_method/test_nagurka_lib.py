from nagurka_fourier_lib import *
import pytest

from math import sin, cos
import pytest_check as check

# N&Y 1988: "This problem can be solved by standard
# linear optimal control methods employing the Hamilton-Jacobi approach via the Riccati equation."

# python -m pytest test_nagurka_lib.py --tb=short -s


def Kirk_A_optimal(t):
    return 7.289*t -6.103 + 6.696*np.exp(-t)-0.593*np.exp(t)

def d_Kirk_A_optimal(t):
    return 7.289 -6.696*np.exp(-t)-0.593*np.exp(t)

def d_d_Kirk_A_optimal(t):
    # not discussed in the paper - just a guess
    return 6.696*np.exp(-t)-0.593*np.exp(t)


def test_Kirk_A_example():
    '''
    Equation 32, Nagurka&Yen 1990, Case A
    '''

    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])

    M = a.shape[0]

    t = np.linspace(epsilon, t_f, 60)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    X = X_(t,P_BCs,a,b,M,t_f)
    d_X = d_X_(t,P_BCs,a,b,M,t_f)
    d_d_X = d_d_X_(t,P_BCs,a,b,M,t_f)

    U = X + d_d_X

    # all of these are perfect except U
    plt.plot(t,X)
    # plt.plot(t,P_(t, p, M))
    # plt.plot(t,
    # plt.plot(t[:-1],np.diff(X)/(1/30))
    plt.plot(t,d_X)
    # plt.show()
    # plt.plot(t[:-2],np.diff(np.diff(X))/((1/30)**2.0))
    plt.plot(t,d_d_X)
    plt.plot(t,Kirk_A_optimal(t))
    plt.plot(t,d_Kirk_A_optimal(t))
    plt.plot(t,d_d_Kirk_A_optimal(t))
    # plt.plot(t,U)
    plt.show()

    # this test fails. the derivatives all work out perfectly, it's just the
    # U control has a completely different shape as fig 1a 1990.
    # must be misunderstanding the U completely somehow.

    #broken
    # assert np.allclose(X,Kirk_A_optimal(t),rtol=0.004)
    # assert np.allclose(d_X,d_Kirk_A_optimal(t),rtol=0.004)
    # assert np.allclose(d_d_X,d_d_Kirk_A_optimal(t),rtol=0.004)

    J = 0.5*integrate.simpson(U**2.0,t)
    assert J == pytest.approx(1.675e1)





#
# def test_Kirk_C_example():
#     '''
#     Equation 42, Nagurka&Yen 1990, Case C
#     '''
#
#     t_f = 2.0
#     a = np.array([4.3833045e-4])
#     b = np.array([1.4015869e-3])
#     M = a.shape[0]
#     t = np.linspace(epsilon, t_f, 60)
#
#     X_t0 = 0.0
#     X_tf = 2.3530569
#     d_X_t0 = 0.0
#     d_X_tf = 2.5293861
#     d_d_X_t0 = 1.3739699
#     d_d_X_tf = 1.9403533
#
#
#
#     p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#
#     X = L_(t,a,b,M,t_f) + P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#     # d_X = d_X_(t,p,a,b,M,t_f)
#     # d_d_X = d_d_X_(t,p,a,b,M,t_f)
#     #
#     U = X + d_d_X
#
#     # plt.plot(t,X)
#     # plt.plot(t,L_(t,a,b,M,t_f))
#
#     # plt.plot(t,P_(t, p, M))
#     # plt.show()
#
#     # assert np.allclose(X,Kirk_A_optimal(t),rtol=0.004)
#     # assert np.allclose(d_X,d_Kirk_A_optimal(t),rtol=0.004)
#     # assert np.allclose(d_d_X,d_d_Kirk_A_optimal(t),rtol=0.004)
#
#     J = 0.5*integrate.simpson(U**2.0,t)
#     assert J == pytest.approx(6.708092)
#
# # def test_transmembrane_compare():
# #
# #     t_f = 2.0
# #     a = np.array([-1.95643e-3])
# #     b = np.array([1.442172e-3])
# #
# #     M = a.shape[0]
# #
# #     t = np.linspace(epsilon, t_f, 60)
# #
# #     virus = default_virus(t)
# #     host_cell = default_host_cell(t)
# #
# #     X_t0 = 0.0
# #     X_tf = 5.0
# #     d_X_t0 = 0.0
# #     d_X_tf = 2.0
# #     d_d_X_t0 = 6.1025137
# #     d_d_X_tf = -3.4798053
# #
# #     p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#

def test_polynomial_BCs():
    t_0 = epsilon
    t_f = 2.0

    a = np.array([1.0, 4.0])
    b = np.array([-0.3,-2])
    M = 2

    t = np.linspace(epsilon, t_f, 100)

    X_t0 = 3.0
    X_tf = 5.0
    d_X_t0 = 4.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    check.almost_equal(P_(epsilon, P_BCs, a, b, M, t_f), X_t0 - L_(epsilon,a,b,M,t_f))
    check.almost_equal(P_(t_f, P_BCs, a, b, M, t_f), X_tf - L_(t_f,a,b,M,t_f))

    check.almost_equal(d_P_(epsilon, P_BCs, a, b, M, t_f), d_X_t0 - d_L_(epsilon,a,b,M,t_f), 1e-4)
    check.almost_equal(d_P_(t_f, P_BCs, a, b, M, t_f), d_X_tf - d_L_(t_f,a,b,M,t_f), 1e-4)

    check.almost_equal(d_d_P_(epsilon, P_BCs, a, b, M, t_f), d_d_X_t0 - d_d_L_(epsilon,a,b,M,t_f), 1e-4)
    check.almost_equal(d_d_P_(t_f, P_BCs, a, b, M, t_f), d_d_X_tf - d_d_L_(t_f,a,b,M,t_f), 1e-4)

#
# def test_L_sin_series():
#     t_f = 1.5
#     M = 1
#
#     # check term a
#     check.almost_equal(1.0, L_(0.0,np.array([1.0]),np.array([0.0]),M,t_f)) # cos 0 is 1
#     check.almost_equal(-1.0, L_(0.5*t_f,np.array([1.0]),np.array([0.0]),M,t_f)) # cos 0 is 1
#     check.almost_equal(1.0, L_(t_f,np.array([1.0]),np.array([0.0]),M,t_f)) # cos 0 is 1
#
#     # d/dx cos 0 is -sin
#     check.almost_equal(0.0, d_L_(0.0,np.array([1.0]),np.array([0.0]),M,t_f))
#     check.almost_equal(-2*pi/t_f, d_L_(0.25*t_f,np.array([1.0]),np.array([0.0]),M,t_f))
#     check.almost_equal(0.0, d_L_(t_f,np.array([1.0]),np.array([0.0]),M,t_f))
#
#     check.almost_equal(4*(pi*pi)/(t_f**2.0), d_d_L_(0.0,np.array([1.0]),np.array([0.0]),M,t_f))
#     check.almost_equal(0.0, d_d_L_(0.25*t_f,np.array([1.0]),np.array([0.0]),M,t_f))
#     check.almost_equal(4*(pi*pi)/(t_f**2.0), d_d_L_(t_f,np.array([1.0]),np.array([0.0]),M,t_f))
#
#     # check term b, sin
#     check.almost_equal(0.0, L_(0.0,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(1.0, L_(0.25*t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(0.0, L_(t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#
#     check.almost_equal(2*pi/t_f, d_L_(0.0,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(0.0, d_L_(0.25*t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(2*pi/t_f, d_L_(t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#
#     check.almost_equal(0.0, d_d_L_(0.0,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(4*(pi*pi)/(t_f**2.0), d_d_L_(0.25*t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#     check.almost_equal(0.0, d_d_L_(t_f,np.array([0.0]),np.array([1.0]),M,t_f))
#
#


'''

Current status:

The time course of U in the first "type A" verification test looks nothing like the plot,
though X is exact

'''
