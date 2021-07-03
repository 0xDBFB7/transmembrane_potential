
'''

Niemann DD. Determination of bang-bang controls for large nonlinear systems n.d.:97.

Uses an endpoint cost by default non-lagrange (or continuous integr

"Initially, one may be tempted to simply sort the switching times
into ascending order before integration. However, the sorting process
destroys the geometric concept of stepping in the negative gradient
direction, and the single variable minimization becomes a rather random process"

#Denote the desired final state vector as Xd.
#A cost function is used to measure the distance between the actual final state X(t) and the Xd
#By determining the partial derivatives of the cost with respect
to the switching times, the switching times can be adjusted to minimize
the state error at the final time.

# 2.11



In Chapter Two, the switching- time vector T and the control vector S were introduced.
Each element of T contains a switching time with the first switch in t1 the second switch in t2, and so on.
The last element of T contains the final time tf. The corresponding
elements of the vector S contain the integer values of the control which
is to be switched. A zero in the final element of S is used to designate the final time.

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



num_switches = 10
T = np.linspace(t_end/10, t_end, num_switches , dtype=np.float128) # switching times vector
# T = np.random.uniform(low=t_end/10, high=t_end, size=(num_switches,)).astype(np.float128)
# S = []


def construct_from_step_functions(T, structure):
    x = np.zeros_like(t)
    x += delta_transmembrane_unit_step(t, structure) # initialize positive
    for idx,T_time in enumerate(T):
        direction = (1.0 - ((idx % 2) * 2.0))
        x -= 2.0 * delta_transmembrane_unit_step(t-T_time, structure) * direction # 1x gets to 0, 2x gets to -1
    return x

def cost_function(T):
    x_v = construct_from_step_functions(T, virus)*1e9
    x_h = construct_from_step_functions(T, host_cell)*1e9
    # return - ((x_v[-1]-(1e-6*1e9))*(x_v[-1]-(1e-6*1e9))) + (x_h[-1]*x_h[-1])
    # return np.sum(x_h*x_h) / np.sum(x_v*x_v)
    return -np.sum(x_v*x_v) + np.sum(x_h*x_h)

# grad_tanh = grad(cost_function)       # Obtain its gradient function
# print(grad_tanh(T))               # Evaluate the gradient at x = 1.0

bounds = [(0, t_end)]*num_switches
from scipy.optimize import minimize,basinhopping

def diagnostics(xk):
    print(xk)

# Tmin = minimize(cost_function, T, method="Powell", options={"disp":True}, callback=diagnostics, bounds=bounds, tol=1e-12).x
tubthumper = basinhopping
minimizer_kwargs = dict(method="Powell", options={"disp":True}, bounds=bounds, callback=diagnostics,  tol=1e-12)
Tmin = tubthumper(cost_function, T, stepsize=t_end/10, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]


plt.plot(t, construct_from_step_functions(Tmin, virus))
plt.plot(t, construct_from_step_functions(Tmin, host_cell))
plt.show()



# plt.draw()
# plt.pause(0.001)
#To ensure that each is greater than the next,


#Now define the (ns+1) x k matrix Ba which
#contains only the columns of B which correspond to active constraints.

# B = np.diag(np.zeros((num_switches,4))) #boundary constraint vector



# input 1(0), X.(0). 11(0), T«-'. S*°>. and 5..,.

# in the paper, they
# Sort T, sort S
# Calculate the cost and the gradient
# Re-sort T and S and the gradient
# however since we only have one control, I don't think we need to do that.
