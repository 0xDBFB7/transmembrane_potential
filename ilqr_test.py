import numpy as np
import theano.tensor as T
import matplotlib.pyplot as plt

from ilqr import iLQR
from ilqr.cost import QRCost
from ilqr.dynamics import AutoDiffDynamics

from transmembrane_lib import *


#https://github.com/anassinator/ilqr/blob/master/examples/rendezvous.ipynb

dt = 0.5e-9  # Discrete time step.
# dt = 0.1
N = 2000  # Number of time steps in trajectory.

t = np.arange(N + 1) * dt

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 50e-9, 5e-9, t)
virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)

a1_v = (virus.a_1)
a2_v = (virus.a_2)
a3_v = (virus.a_3)
b1_v = (virus.b_1)
b2_v = (virus.b_2)
b3_v = (virus.b_3)
R_v = (virus.R)


a1_h = (host_cell.a_1)
a2_h = (host_cell.a_2)
a3_h = (host_cell.a_3)
b1_h = (host_cell.b_1)
b2_h = (host_cell.b_2)
b3_h = (host_cell.b_3)
R_h = (host_cell.R)



def on_iteration(iteration_count, xs, us, J_opt, accepted, converged):
    J_hist.append(J_opt)
    info = "converged" if converged else ("accepted" if accepted else "failed")
    print("iteration", iteration_count, info, J_opt)

xEq_v = T.dscalar("xEq_v")

x0_v = T.dscalar("x0_v")
x1_v = T.dscalar("x1_v")


x0_h = T.dscalar("x0_h")
x1_h = T.dscalar("x1_h")

u0 = T.dscalar("u0")
u1 = T.dscalar("u1")
u2 = T.dscalar("u2")

x_inputs = [
    xEq_v,
    x0_v,
    x1_v,
    x0_h,
    x1_h,
    u0,
    u1
]

u_inputs = [
    u2
]


# Discrete dynamics model definition.
f = T.stack([
    xEq_v + (x0_v - 1.0),
    x0_v + x1_v * dt,
    x1_v + ((R_v*a1_v*u2 + R_v*a2_v*u1 + R_v*a3_v*u0 - b2_v*x1_v - b3_v*(x0_v))/b1_v) * dt,
    # x2_v + ,

    x0_h + (x1_h * dt),
    x1_h + ((R_h*a1_h*u2 + R_h*a2_h*u1 + R_h*a3_h*u0 - b2_h*x1_h - b3_h*x0_h)/b1_h) * dt,
    # x2_h + ,

    u0 + u1 * dt,
    u1 + u2 * dt
])


dynamics = AutoDiffDynamics(f, x_inputs, u_inputs)



# Q = np.eye(dynamics.state_size)#state error

# cost = transpose(x) * Q * x + transpose(u) * R * u


Q = np.zeros((dynamics.state_size,dynamics.state_size))

Q[0,0] = 1 # xExp^2
Q[3,3] = 1 # xExp^2

# Q[5,5] = 0.
#One of the essential features of LQR is that Q should be a
#symmetric positive semidefinite matrix and R should be a positive definite matrix.
# Q = 1
# Q[0, 0] = Q[0, 0] = -1
# Q[0, 0] = 1



R = 0.0 * np.eye(dynamics.action_size)#control error

cost = QRCost(Q, R)

x0 = np.array([0, 0, 0, 0, 0, 0, 0])  # Initial state.

# Random initial action path.
us_init = np.random.uniform(-1, 1, (N, dynamics.action_size))

J_hist = []
ilqr = iLQR(dynamics, cost, N)
xs, us = ilqr.fit(x0, us_init, on_iteration=on_iteration)

# test run
#
def test_run():
    us_init = np.exp(-(((t[0:-1]-((N//4)*dt))**2.0)/(2.0*((1e-9)**2.0)))) # for simulation
    us_init = np.diff(us_init, prepend=0.0)
    us_init = np.diff(us_init,  prepend=0.0)


    (xs, F_x, F_u, L, L_x, L_u, L_xx, L_ux, L_uu, F_xx, F_ux,
     F_uu) = ilqr._forward_rollout(x0, us_init)




x_v = xs[:, 1]
x_h = xs[:, 3]

u = xs[:, 5]


plt.subplot(1,3,1)

plt.plot(t, x_v, "r")
plt.subplot(1,3,2)
plt.plot(t, x_h, "b")
# _ = plt.plot(t, u0, "g")

plt.subplot(1,3,3)
plt.plot(t, u, "r")




plt.show()
plt.figure()
ideal_values = u
plt.plot(t, convolve_output(ideal_values, host_cell, dt))
plt.plot(t, convolve_output(ideal_values, virus, dt))

print((np.max(convolve_output(ideal_values, host_cell, dt) * 1e6) - np.min(convolve_output(ideal_values, host_cell, dt) * 1e6))
            /(np.max(convolve_output(ideal_values, virus, dt) * 1e6) - np.min(convolve_output(ideal_values, virus, dt) * 1e6)))


plt.show()
