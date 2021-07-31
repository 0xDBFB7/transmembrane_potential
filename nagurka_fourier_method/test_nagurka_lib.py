from naguka_fourier_lib import *
import pytest


def Kirk_A_optimal(t):
    return 7.289*t -6.103 + 6.696*np.exp(-t)-0.593*np.exp(t)

def d_Kirk_A_optimal(t):
    return 7.289 -6.696*np.exp(-t)-0.593*np.exp(t)

def test_Kirk_example():
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
    d_d_X_t0 = 6.102
    d_d_X_tf = -3.479


    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    X = P_(t,p,M) + L_(t,a,b,M,t_f)
    d_X = d_P_(t,p,M) + d_L_(t,a,b,M,t_f)
    d_d_X = d_d_P_(t,p,M) + d_d_L_(t,a,b,M,t_f)



    print(P_(epsilon,p,M) + L_(epsilon,a,b,M,t_f))
    print(P_(2.0,p,M) + L_(2.0,a,b,M,t_f))

    print(d_P_(epsilon,p,M) + d_L_(epsilon,a,b,M,t_f))
    print(d_P_(2.0,p,M) + d_L_(2.0,a,b,M,t_f))

    print(d_d_P_(epsilon,p,M) + d_d_L_(epsilon,a,b,M,t_f))
    print(d_d_P_(2.0,p,M) + d_d_L_(2.0,a,b,M,t_f))
    # The performance index was evaluated by means of Simpson's composite
    # integral formula with a step size of 1/30 (consistent unit of
    # time).

    # all of these are perfect except U, need asserts 
    # plt.plot(t,U)
    # plt.plot(t,X)
    # plt.plot(t[:-1],np.diff(X)/(1/30))
    # plt.plot(t,d_X)
    # plt.show()
    # plt.plot(t[:-2],np.diff(np.diff(X))/((1/30)**2.0))
    # plt.plot(t,d_d_X)
    # plt.plot(t,Kirk_A_optimal(t))
    # plt.plot(t,d_Kirk_A_optimal(t))
    # plt.show()


    U = X + d_d_X

    # this test fails. the derivatives all seem to work out, it's just the
    # U control has a completely different shape!
    # must be misunderstanding the U

    J = 0.5*integrate.simpson(U**2.0,t)
    assert J == pytest.approx(1.675e1)
