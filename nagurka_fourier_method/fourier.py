
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

from autograd import grad    # The only autograd function you may ever need
import sys
sys.path.append('../')
from transmembrane_lib import *
import matplotlib.pyplot as plt


# actually, we don't even need to integrate - just a sum of step functions will do!

t_end = 5e-10
num_time_nodes = 5000
t = np.linspace(0, t_end, num_time_nodes, dtype=np.float128) # integration span time vector

virus = default_virus(t)
host_cell = default_host_cell(t)





# Tmin = minimize(cost_function, T, method="Powell", options={"disp":True}, callback=diagnostics, bounds=bounds, tol=1e-12).x
tubthumper = basinhopping
minimizer_kwargs = dict(method="Powell", options={"disp":True}, bounds=bounds, callback=diagnostics,  tol=1e-12)
Tmin = tubthumper(cost_function, T, stepsize=t_end/10, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]
