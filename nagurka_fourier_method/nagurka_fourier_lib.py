
'''

Fourier-Based Optimal Control of Nonlinear Dynamic Systems, 1990 -
Journal of Dynamic Systems, Measurement, and Control,


F[..] = U(t)

F is the set of ordinary differential equations used to describe
the system behavior

Here, E represents the np.cost associated with the terminal states and G represents the np.cost associated with the trajectory


The optimal trajectory, X*(/), X*(0, and X*0, is defined
as the admissible trajectory that minimizes the performance
index, J.


"np.since the derivatives of the
generalized coordinates are obtained by direct analytical differention of equation (13), the system equations (1) are treated
as algebraic equations in evaluating the control variables. The
computational scheme of the proposed approach is therefore an inverse dynamic method.
As a result, no integration of differential equations (such as state and np.costate equations) is required.
The computational np.cost is thus significantly reduced."



x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))


#Note: N&Y 1988 vs 1990 uses confunp.singly different notation. beware!

There are a number of different versions of the papers by N&Y.


'''


# from autograd.scipy.integrate import odeint
from scipy.integrate import odeint
# import autograd.numpy as np  # Thinly-wrapped numpy
import numpy as np
import scipy.integrate as integrate
from icecream import ic
from math import pi
#maybe jax someday?




import matplotlib.pyplot as plt
epsilon = 1e-15



# t_ = 5e-10
# num_time_nodes = 5000
# t = np.linspace(0, t_, num_time_nodes, dtype=np.float128) # integration span time vector


# The problem is solved X_v = P + lambda (that is, we solve for one time-course of the virus,
# then
# Xv -> U (obtain the control required to produce that time course
# then U ->
# this is a little clunky, not sure why I chose this, but it's
# along the lines of what's recommed by N&Y.

# M = a.shape[0]+1 # is this right?


def L_(t, a, b, M, t_f):
    L = t*0.0
    m = np.arange(1, M+1)
    for i in m:
        v = i*(2*pi / t_f)
        L += a[i-1]*np.cos(v*t) + b[i-1]*np.sin(v*t)
    return L


# d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, M, t_f), t)
# d_d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, M, t_f), t)

def d_L_(t, a, b, M, t_f):
    L = t*0.0
    m = np.arange(1, M+1)
    for i in m:
        v = i*(2*pi / t_f)
        L += v * -1.0 * a[i-1]*np.sin(v*t) + v * b[i-1]*np.cos(v*t)
    return L


def d_d_L_(t, a, b, M, t_f):
    L = t*0.0
    m = np.arange(1, M+1)
    for i in m:
        v = i*(2*pi / t_f)
        L += (v**2) * -1.0 * a[i-1]*np.cos(v*t) + (v**2) * -1.0 * b[i-1]*np.sin(2*pi*(m*t)/t_f)
    return L


def X_(t,p,a,b,M,t_f):
    return P_(t, p, a, b, M, t_f) + L_(t, a, b, M, t_f)


def d_X_(t,p,a,b,M,t_f):
    return d_P_(t, p, a, b, M, t_f) + d_L_(t, a, b, M, t_f)


def d_d_X_(t,p,a,b,M,t_f):
    return d_d_P_(t, p, a, b, M, t_f) + d_d_L_(t, a, b, M, t_f)


# should the fourier series be always be *zero* at the edges, or simply *the same*?

def X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
    """
    Moves polynomial boundary conditions so that the fourier series is
    always zero at the edges for convergence
    """
    return (X_t0 - L_(epsilon, a, b, M, t_f),
            d_X_t0 - d_L_(epsilon, a, b, M, t_f),
            d_d_X_t0 - d_d_L_(epsilon, a, b, M, t_f),

            X_tf - L_(t_f, a, b, M, t_f),
            d_X_tf - d_L_(t_f, a, b, M, t_f),
            d_d_X_tf - d_d_L_(t_f, a, b, M, t_f))


def P_BCs_to_p_coefficients(P_BCs, t_f):
    """
    Converts the polynomial boundary conditions that must be enforced into
    the polynomial coefficients
    """
    P_t0, d_P_t0, d_d_P_t0, P_tf, d_P_tf, d_d_P_tf = P_BCs

    p_0 = P_t0
    p_1 = d_P_t0
    p_2 = d_d_P_t0/2
    p_3 = (t_f**2*(-3*d_d_P_t0 + d_d_P_tf) - 4*t_f*(3*d_P_t0 + 2*d_P_tf) - 20*P_t0 + 20*P_tf)/(2*t_f**3)
    p_4 = (t_f**2*(3*d_d_P_t0 - 2*d_d_P_tf)/2 + t_f*(8*d_P_t0 + 7*d_P_tf) + 15*P_t0 - 15*P_tf)/t_f**4
    p_5 = (t_f**2*(-d_d_P_t0 + d_d_P_tf) - 6*t_f*(d_P_t0 + d_P_tf) - 12*P_t0 + 12*P_tf)/(2*t_f**5)

    return (p_0, p_1, p_2, p_3, p_4, p_5)


def P_(t, P_BCs, a, b, M, t_f):
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = P_BCs_to_p_coefficients(P_BCs, t_f)

    return p_0 + p_1*t + p_2*t**2 + p_3*t**3 + p_4*t**4 + p_5*t**5


# d_P_(t, P_BCs, a, b, M, t_f) = ForwardDiff.derivative(n -> P_(n, P_BCs, a, b, M, t_f), t)
# d_d_P_(t, P_BCs, a, b, M, t_f) = ForwardDiff.derivative(n -> d_P_(n, P_BCs, a, b, M, t_f), t)

def d_P_(t, P_BCs, a, b, M, t_f):
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = P_BCs_to_p_coefficients(P_BCs, t_f)

    return p_1 + 2*p_2*t + 3*p_3*t**2 + 4*p_4*t**3 + 5*p_5*t**4


def d_d_P_(t, P_BCs, a, b, M, t_f):
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = P_BCs_to_p_coefficients(P_BCs, t_f)

    return 2*p_2 + 6*p_3*t + 12*p_4*t**2 + 20*p_5*t**3
