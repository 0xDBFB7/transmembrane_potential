
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
import autograd.numpy as np  # Thinly-wrapped numpy
#maybe jax someday?

from autograd import grad    # The only autograd function you may ever need
import sys
sys.path.append('../')
from transmembrane_lib import *


# actually, we don't even need to integrate - just a sum of step functions will do!

t_end = 10e-6
t = np.linspace(0, t_end, 1000) # integration span time vector

virus = default_virus(t)
host_cell = default_host_cell(t)



num_switches = 10
T = np.linspace(0, t_end, num_switches) # switching times vector
# S = []


def cost_function(T):
    x_v = np.zeros_like(T)
    x_h = np.zeros_like(T)
    for idx,T_time in enumerate(T):
        direction = (1.0 - ((idx % 2) * 2.0))
        x_v += delta_transmembrane_unit_step(t-T_time, virus) * direction
        x_h += delta_transmembrane_unit_step(t-T_time, host_cell) * direction

    return -(x_v*x_v) + (x_h*x_h)

grad_tanh = grad(cost_function)       # Obtain its gradient function
grad_tanh(T)               # Evaluate the gradient at x = 1.0


#To ensure that each is greater than the next,


#Now define the (ns+1) x k matrix Ba which
#contains only the columns of B which correspond to active constraints.

B = np.diag(np.zeros(num_switches,4)) #boundary constraint vector


# input 1(0), X.(0). 11(0), T«-'. S*°>. and 5..,.

# in the paper, they
# Sort T, sort S
# Calculate the cost and the gradient
# Re-sort T and S and the gradient
# however since we only have one control, I don't think we need to do that.
