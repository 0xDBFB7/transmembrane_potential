from nagurka_fourier_lib import *
import pytest

from nagurka_fourier_lib_restate import P_restate, P_restate_coeffs
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


def test_Kirk_A_example_statement_one():
    '''
    Equation 32, Nagurka&Yen 1990, Case A
    '''

    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])

    M = a.shape[0]

    t = np.linspace(epsilon, t_f, 60)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    X = X_(t,p,a,b,M,t_f)
    d_X = d_X_(t,p,a,b,M,t_f)
    d_d_X = d_d_X_(t,p,a,b,M,t_f)

    U = X + d_d_X

    # all of these are perfect except U, need asserts
    # plt.plot(t,X)
    # plt.plot(t,P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M))
    # plt.plot(t,
    # plt.plot(t[:-1],np.diff(X)/(1/30))
    # plt.plot(t,d_X)
    # plt.show()
    # plt.plot(t[:-2],np.diff(np.diff(X))/((1/30)**2.0))
    # plt.plot(t,d_d_X)
    # plt.plot(t,Kirk_A_optimal(t))
    # # plt.plot(t,d_Kirk_A_optimal(t))
    # # plt.plot(t,d_d_Kirk_A_optimal(t))
    # plt.plot(t,U)
    # plt.show()

    # this test fails. the derivatives all work out perfectly, it's just the
    # U control has a completely different shape as fig 1a 1990.
    # must be misunderstanding the U completely somehow.

    print(np.max(np.abs(X-Kirk_A_optimal(t))))
    assert np.allclose(X,Kirk_A_optimal(t),rtol=0.004)
    assert np.allclose(d_X,d_Kirk_A_optimal(t),rtol=0.004)
    assert np.allclose(d_d_X,d_d_Kirk_A_optimal(t),rtol=0.004)

    J = 0.5*integrate.simpson(U**2.0,t)
    assert J == pytest.approx(1.675e1)

def test_transmembrane_compare():

    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])

    M = a.shape[0]

    t = np.linspace(epsilon, t_f, 60)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

#
def test_BCs():

    t_f = 1.0

    a = np.array([1.0, 1.0])
    b = np.array([1.0, -1.0])

    M = a.shape[0]

    t = np.linspace(0, t_f, 100)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    X_t0 = 1.0 # is it possible that this is overconstrained?
    X_tf = 2.0
    d_X_t0 = 0.0
    d_X_tf = 0.0
    d_d_X_t0 = 0.0
    d_d_X_tf = 0.0

    # Here X, (2) is selected as being free while ^,(2) is calculated according to Xl(2) = [15 - X(2)]/5

    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    # plt.plot(t,P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M))
    # plt.plot(t,P_(t,p,M))
    # plt.plot(t,L_(t,a,b,M,t_f))
    # plt.plot(t,L_(t,a,b,M,t_f) + P_(t,p,M))

    # plt.plot(t, P_restate(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M) - P_(t,p,M))
    # plt.show()

    # assert X_(epsilon,p,a,b,M,t_f) == pytest.approx(X_t0)
    # assert X_(t_f,p,a,b,M,t_f) == pytest.approx(X_tf)

    assert P_restate(t_f, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M) +\
                        L_(t_f,a,b,M,t_f) == pytest.approx(X_tf)


def test_P_coeffs():
    t_f = 2.0

    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])
    a_ = -1.95643e-3
    b_ = 1.442172e-3
    M = 1

    t = np.linspace(epsilon, t_f, 100)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    check.almost_equal(P_(0.0,p,M), X_t0 - L_(0.0,a,b,M,t_f))
    check.almost_equal(P_(t_f,p,M), X_tf - L_(t_f,a,b,M,t_f))

    check.almost_equal(d_P_(0.0,p,M), d_X_t0 - d_L_(0.0,a,b,M,t_f), 1e-4)
    check.almost_equal(d_P_(t_f,p,M), d_X_tf - d_L_(t_f,a,b,M,t_f), 1e-4)

    check.almost_equal(d_d_P_(0.0,p,M), d_X_t0 - d_d_L_(0.0,a,b,M,t_f), 1e-4)
    check.almost_equal(d_d_P_(t_f,p,M), d_X_tf - d_d_L_(t_f,a,b,M,t_f), 1e-4)


def test_P_restate_coeffs():
    t_f = 2.0

    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])
    a_ = -1.95643e-3
    b_ = 1.442172e-3
    M = 1

    t = np.linspace(epsilon, t_f, 100)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    # X_t0 = 1.0 # is it possible that this is overconstrained?
    # X_tf = 1.0
    # d_X_t0 = 0
    # d_X_tf = 0
    # d_d_X_t0 = 0
    # d_d_X_tf = 0

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = P_restate_coeffs(t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    # https://gamma.sympy.org/
    # diff(cos(2*pi*t/t_f), t)
    # d_L_t0 = ((a_ * -2 * pi * sin(2 * pi * 0.0 / t_f)/t_f) + (b_ * 2 * pi * cos(2 * pi * 0.0 / t_f)/t_f))
    # d_L_tf = ((a_ * -2 * pi * sin(2 * pi * t_f / t_f)/t_f) + (b_ * 2 * pi * cos(2 * pi * t_f / t_f)/t_f))

    check.almost_equal(p0, X_t0 - L_(0.0,a,b,M,t_f))
    check.almost_equal(P_restate(0.0, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M), X_t0 - L_(0.0,a,b,M,t_f))

    check.almost_equal(pf, X_tf - L_(t_f,a,b,M,t_f))
    check.almost_equal(P_restate(t_f, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M), X_tf - L_(t_f,a,b,M,t_f))

    check.almost_equal(d_p0, d_X_t0 - d_L_(0.0,a,b,M,t_f), 1e-4)
    check.almost_equal(d_pf, d_X_tf - d_L_(t_f,a,b,M,t_f), 1e-4)
    check.almost_equal(d_d_p0, d_X_t0 - d_d_L_(0.0,a,b,M,t_f), 1e-4)
    check.almost_equal(d_d_pf, d_X_tf - d_d_L_(t_f,a,b,M,t_f), 1e-4)


'''

Current status:

Something to do with the b term seems to be broken somehow.
The time course of U in the first "type A" verification test looks nothing like the plot,
though X is exact

'''
