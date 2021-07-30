from naguka_fourier_lib import *
import pytest



def test_Kirk_example():



    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])

    M = a.shape[0]

    t = np.linspace(0, t_f, 100)

    virus = default_virus(t)
    host_cell = default_host_cell(t)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 0.0
    d_d_X_tf = 0.0



    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    X = P_(t,p,M) + L_(t,a,b,M,t_f)
    d_d_X = d_d_P_(t,p,M) + d_d_L_(t,a,b,M,t_f)

    # The performance index was evaluated by means of Simpson's composite
    # integral formula with a step size of 1/30 (consistent unit of
    # time).

    U = X + d_d_X
    J = integrate.simpson(U**2.0,t)
    print(J)
    assert J == pytest.approx(1.675e1)
