
'''

Fourier-Based Optimal Control of Nonlinear Dynamic Systems, 1990 -
Journal of Dynamic Systems, Measurement, and Control,


F[..] = U(t)

F is the set of ordinary differential equations used to describe
the system behavior

Here, E represents the cost associated with the terminal states and G represents the cost associated with the trajectory


The optimal trajectory, X*(/), X*(0, and X*0, is defined
as the admissible trajectory that minimizes the performance
index, J.


"Since the derivatives of the
generalized coordinates are obtained by direct analytical differention of equation (13), the system equations (1) are treated
as algebraic equations in evaluating the control variables. The
computational scheme of the proposed approach is therefore an inverse dynamic method.
As a result, no integration of differential equations (such as state and costate equations) is required.
The computational cost is thus significantly reduced."



x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))


#Note: N&Y 1988 vs 1990 uses confusingly different notation. beware!

There are a number of different versions of the papers by N&Y.


'''


# from autograd.scipy.integrate import odeint
from scipy.integrate import odeint
# import autograd.numpy as np  # Thinly-wrapped numpy
import numpy as np
import scipy.integrate as integrate
from icecream import ic
#maybe jax someday?

# from autograd import grad    # The only autograd function you may ever need

import sys
sys.path.append('../')
from transmembrane_lib import *

import matplotlib.pyplot as plt
epsilon = 1e-20



# t_end = 5e-10
# num_time_nodes = 5000
# t = np.linspace(0, t_end, num_time_nodes, dtype=np.float128) # integration span time vector


# The problem is solved X_v = P + lambda (that is, we solve for one time-course of the virus,
# then
# Xv -> U (obtain the control required to produce that time course
# then U ->
# this is a little clunky, not sure why I chose this, but it's
# along the lines of what's recommended by N&Y.

# M = a.shape[0]+1 # is this right?
#
# class parametrization():

# def __init__(s, t, X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):


def P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M):
    '''
    Eqs 22 to 27, Nagurka&Yen 1990
    1988 has an alternate setup
    '''

    m = np.arange(1, M+1)
    p = np.zeros((6))

    p[0] = X_t0 - np.sum(a)
    #
    p[1] = d_X_t0 - (2 * pi / t_f) * np.sum(m * b)
    #
    p[2] = (d_d_X_t0/2.0) + (2 * (pi**2.0) / (t_f**2.0)) * np.sum((m**2.0) * a)
    #
    p[3] = ((10.0*(X_tf - X_t0)) + (20.0*pi*np.sum(m * b)) - 4.0*(pi**2.0)*np.sum((m**2.0)*a))/(t_f**3.0)
    p[3] += -(6*d_X_t0 + 4*d_X_tf)/(t_f**2.0)
    p[3] += -((3.0/2.0)*d_d_X_t0 - (0.5*d_d_X_tf))/(t_f)
    #
    p[4] = (15*(X_t0 - X_tf) - (30*pi*np.sum(m*b))+(2*(pi**2.0)*np.sum((m**2.0)*a)))/(t_f**4.0)
    p[4] += (8.0*d_X_t0 + 7*d_X_tf)/(t_f**3.0)
    p[4] += (3/2.0*d_d_X_t0 - d_d_X_tf)/(t_f**2.0)
    #
    p[5] = (6.0*(X_tf-X_t0) + (4.0*pi*np.sum(m*b)))/(t_f**5.0)
    p[5] += -(3*d_X_t0 + 3*d_X_tf)/(t_f**4.0) - 0.5*(d_d_X_t0 - d_d_X_tf)/(t_f**3.0)

    return p

def P_(t, p, M):
    P = t*0.0 # supports both float and array type
    for j in range(0,6):
        P += p[j]*(t**j)
    return P

def d_P_(t, p, M):
    P = t*0.0
    for j in range(0,6):
        P += (j*(t**j)*p[j])

    if(isinstance(t,np.ndarray)):
        P = np.divide(P, t, out=np.zeros_like(P), where=(t)!=0)
    else:
        P = P / t

    return P

def d_d_P_(t, p, M):
    P = t*0.0

    for j in range(0,6):
        P += j*(t**j)*(j - 1)*p[j]

    if(isinstance(t,np.ndarray)):
        P = np.divide(P, (t**2), out=np.zeros_like(P), where=(t**2)!=0)
    else:
        P = P / (t**2.0)


    return P

    # if(isinstance(t,np.ndarray)):

def L_(t,a,b,M,t_f): # lambda
    L = t*0.0
    for m in range(1,M+1):
        L+= np.sin(2*pi*m*t/t_f)*b[m-1]
    for m in range(1,M+1):
        L+= np.cos(2*pi*m*t/t_f)*a[m-1]
    return L

def d_L_(t,a,b,M,t_f):
    L = t*0.0
    for m in range(1,M+1):
        L += -2*pi*m*np.sin(2*pi*m*t/t_f)*a[m-1]/t_f
    for m in range(1,M+1):
        L += 2*pi*m*np.cos(2*pi*m*t/t_f)*b[m-1]/t_f
    return L

#sympy lambdify would probably be more straightforward
def d_d_L_(t,a,b,M,t_f):
    L = t*0.0

    for m in range(1,M+1):
        L+= m**2*np.sin(2*pi*m*t/t_f)*b[m-1] * -4*(pi**2.0)/(t_f**2.0)
    for m in range(1,M+1):
        L+= m**2*np.cos(2*pi*m*t/t_f)*a[m-1] * -4*(pi**2.0)/(t_f**2.0)

    return L

def X_(t,p,a,b,M,t_f):
    return P_(t,p,M) + L_(t,a,b,M,t_f)

def d_X_(t,p,a,b,M,t_f):
    return d_P_(t,p,M) + d_L_(t,a,b,M,t_f)

def d_d_X_(t,p,a,b,M,t_f):
    return d_d_P_(t,p,M) + d_d_L_(t,a,b,M,t_f)



# def Xv_to_
