import numpy as np
import theano.tensor as T
import matplotlib.pyplot as plt

from ilqr import iLQR
from ilqr.cost import QRCost
from ilqr.dynamics import AutoDiffDynamics
from ilqr.dynamics import FiniteDiffDynamics

from ilqr.dynamics import BatchAutoDiffDynamics

from transmembrane_lib import *


#https://github.com/anassinator/ilqr/blob/master/examples/rendezvous.ipynb

T0 = 1e-8
U0 = 1.0
X0 = 1e-9

N = 3000  # Number of time steps in trajectory.
end = 1e-5 / T0
dt = end/N  # Discrete time step.

# dt = 0.1



t = np.arange(N + 1) * dt

virus = default_virus(t)
host_cell = default_host_cell(t)

def on_iteration(iteration_count, xs, us, J_opt, accepted, converged):
    J_hist.append(J_opt)
    info = "converged" if converged else ("accepted" if accepted else "failed")
    print("iteration", iteration_count, info, J_opt)


x0_v = T.dscalar("x0_v")
x1_v = T.dscalar("x1_v")


x0_h = T.dscalar("x0_h")
x1_h = T.dscalar("x1_h")

u0 = T.dscalar("u0")
u1 = T.dscalar("u1")
u2 = T.dscalar("u2")

x_inputs = [
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
    N + d_pore_density(x0_v, N, pore_N0, pore_alpha, pore_q, pore_V_ep)
    x0_v + (x1_v * dt),
    x1_v + (virus.alpha*u2 + virus.beta*u1 + virus.gamma*u0 - virus.phi*x1_v - virus.xi*x0_v)) * dt,
    # x2_v + ,

    N + d_pore_density(x0_v, N, pore_N0, pore_alpha, pore_q, pore_V_ep)
    x0_h + (x1_h * dt),
    x1_h + (((U0 / (T0**2))*alpha_h*u2 + (U0 / T0)*beta_h*u1 + gamma_h*U0*u0 - phi_h*(X0 / T0)*x1_h - xi_h*X0*x0_h)/(X0 / (T0**2))) * dt,
    # x2_h + ,

    u0 + (u1 * dt),
    u1 + (u2 * dt)
])


dynamics = AutoDiffDynamics(f, x_inputs, u_inputs)
# dynamics = FiniteDiffDynamics(f, 6, 1)
# dynamics = BatchAutoDiffDynamics(f, state_size, action_size)

# Q = np.eye(dynamics.state_size)#state error

# cost = transpose(x) * Q * x + transpose(u) * R * u


Q = np.zeros((dynamics.state_size,dynamics.state_size))

Q[0,0] = 1 # xExp^2
Q[2,2] = 1 # xExp^2

# Q[4,4] = 0.1
#One of the essential features of LQR is that Q should be a
#symmetric positive semidefinite matrix and R should be a positive definite matrix.
# Q = 1
# Q[0, 0] = Q[0, 0] = -1
# Q[0, 0] = 1



R = 0.0 * np.eye(dynamics.action_size)#control error

x0 = np.array([0, 0, 0, 0, 0, 0])  # Initial state.
x_goal = np.array([1, 0, 0, 0, 0, 0])  # Initial state.

cost = QRCost(Q, R, x_goal=x_goal)



# Random initial action path.
us_init = np.random.uniform(-1, 1, (N, dynamics.action_size))

J_hist = []
ilqr = iLQR(dynamics, cost, N)
xs, us = ilqr.fit(x0, us_init, n_iterations=200, on_iteration=on_iteration)

# test run
#
def test_run():
    us_init = np.exp(-(((t[0:-1]-((N//4)*dt))**2.0)/(2.0*((1e-9)**2.0)))) # for simulation
    us_init = np.diff(us_init, prepend=0.0)
    us_init = np.diff(us_init,  prepend=0.0)


    (xs, F_x, F_u, L, L_x, L_u, L_xx, L_ux, L_uu, F_xx, F_ux,
     F_uu) = ilqr._forward_rollout(x0, us_init)




x_v = xs[:, 0]
x_h = xs[:, 2]

u = xs[:, 4]


plt.subplot(1,3,1)

plt.plot(t, x_v, "r")
plt.plot(t, xs[:, 0], "r")

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

print((np.max(convolve_output(ideal_values, host_cell, dt) ) - np.min(convolve_output(ideal_values, host_cell, dt)))
            /(np.max(convolve_output(ideal_values, virus, dt) ) - np.min(convolve_output(ideal_values, virus, dt))))

print( np.abs(np.sum(convolve_output(ideal_values, virus, dt))) / np.abs(np.sum(convolve_output(ideal_values, host_cell, dt))))

plt.show()
