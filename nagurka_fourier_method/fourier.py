
'''


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

'''


from autograd.scipy.integrate import odeint
# import autograd.numpy as np  # Thinly-wrapped numpy
import numpy as np
#maybe jax someday?

# from autograd import grad    # The only autograd function you may ever need

import sys
sys.path.append('../')
from transmembrane_lib import *

import matplotlib.pyplot as plt



virus = default_virus(t)
host_cell = default_host_cell(t)


# t_end = 5e-10
# num_time_nodes = 5000
# t = np.linspace(0, t_end, num_time_nodes, dtype=np.float128) # integration span time vector


# The problem is solved X_v = P + lambda (that is, we solve for one time-course of the virus,
# then
# Xv -> U (obtain the control required to produce that time course
# then U ->
# this is a little clunky, not sure why I chose this, but it's
# along the lines of what's recommended by N&Y.

def P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b):
    '''
    Eqs 22 to 27, Nagurka&Yen
    '''
    M = a.shape[0]+1 # is this right?
    m = np.arange(1, M)
    p = np.zeros((5))
    p[0] = X_t0 - np.sum(a)
    p[1] = d_X_t0 - (2 * pi / t_f) * np.sum(m * b)
    p[2] = (d_d_X_t0/2.0) - (2 * (pi**2.0) / (t_f**2.0)) * np.sum((m**2.0) * a)

    p[3] = ((10.0*(X_tf - X_t0)) + (20.0*pi*np.sum(m * b)) - 4.0*(pi**2.0)*np.sum((m**2.0)*a))/(t_f**3.0)
    p[3] += -(6*d_X_t0 + 4*d_X_tf)/(t_f**2.0)
    p[3] += -(((3.0/2.0)*d_X_t0 - (0.5*d_d_X_tf))/(t_f)

    p[4] = (15*(X_t0 - X_tf) - (30*pi*np.sum(m*b))+(2*(pi**2.0)*np.sum((m**2.0)*a)))/(t_f**4.0)
    p[4] += (8.0*d_X_t0 + 7*d_X_tf)/(t_f**3.0)
    p[4] += (3/2.0*d_d_X_t0 - d_d_X_tf)/(t_f**2.0)

    p[5] = (6.0*(X_tf-X_t0) + (4.0*pi*np.sum(m*b)))/(t_f**5.0)
    p[5] = -(3*d_X_t0 + 3*d_X_tf)/(t_f**4.0) - 0.5*(d_d_X_t0 - d_d_X_tf)/(t_f**3.0)
    
    return p
# unit test - take numerical derivative and compare?

def P():


Sum(t**j*p[j], (j, 0, 5)) + Sum(sin(2*pi*m*t/t_f)*b[m], (m, 1, M)) + Sum(cos(2*pi*m*t/t_f)*a[m], (m, 1, M))
Sum(j*t**j*p[j]/t, (j, 0, 5)) + Sum(-2*pi*m*sin(2*pi*m*t/t_f)*a[m]/t_f, (m, 1, M)) + Sum(2*pi*m*cos(2*pi*m*t/t_f)*b[m]/t_f, (m, 1, M))
-4*pi**2*Sum(m**2*sin(2*pi*m*t/t_f)*b[m], (m, 1, M))/t_f**2 - 4*pi**2*Sum(m**2*cos(2*pi*m*t/t_f)*a[m], (m, 1, M))/t_f**2 + Sum(j*t**j*(j - 1)*p[j], (j, 0, 5))/t**2



# Tmin = minimize(cost_function, T, method="Powell", options={"disp":True}, callback=diagnostics, bounds=bounds, tol=1e-12).x
# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Powell", options={"disp":True}, bounds=bounds, callback=diagnostics,  tol=1e-12)
# Tmin = tubthumper(cost_function, T, stepsize=t_end/10, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]
