from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os
from scipy import interpolate
from gekko import GEKKO

def normalized_gaussian_pulse(t,fwhm):
    sigma = fwhm/2.355
    return np.exp(-((t**2.0)/(2.0*(sigma**2.0))))


T0 = 1e-8


sim_t = np.loadtxt( 'PSOPT/build/t.dat', dtype=np.float128) * T0
sim_t = sim_t[:-1] # duplicate!

t0 = 0
tstop = sim_t[-1]

dt = 0.00005e-9
t = np.linspace(t0, tstop, int(tstop/dt), dtype=np.float128)

# dt = (t[1] - t[0])

sim_values = np.loadtxt( 'PSOPT/build/u0.dat', dtype=np.float128 )[:-1]


# flinear = interpolate.interp1d(x, sim_values)
# fcubic = interpolate.interp1d(sim_t, sim_values, kind="quadratic")
interpolant = interpolate.InterpolatedUnivariateSpline(sim_t, sim_values, k=5)

control_input = interpolant(np.float64(t))
host_cell = default_host_cell(t)
virus = default_virus(t)


m = GEKKO() # initialize gekko
m.options.MAX_ITER = 200
nt = 800
U0 = 1.0
X0 = 1e-6
end = tstop
m.time = np.linspace(0,end,nt)

gekko_input = np.float64(interpolant(np.float64(m.time)))
u0 = m.Param(value=gekko_input)
u1_ = np.diff(gekko_input, prepend=gekko_input[0])/((end/nt)/T0)
u1 = m.Param(value=u1_)
u2_ = np.diff(u1_, prepend=u1_[0])/((end/nt)/T0)
u2 = m.Param(value=u2_)

x0_v = m.Var(value=0)
x1_v = m.Var(value=0)
x2_v = m.Var(value=0)

x0_h = m.Var(value=0)
x1_h = m.Var(value=0)
x2_h = m.Var(value=0)

# t = m.Param(value=m.time)

alpha_h = m.Const(host_cell.alpha)
beta_h = m.Const(host_cell.beta)
gamma_h = m.Const(host_cell.gamma)
phi_h = m.Const(host_cell.phi)
xi_h = m.Const(host_cell.xi)

alpha_v = m.Const(virus.alpha)
beta_v = m.Const(virus.beta)
gamma_v = m.Const(virus.gamma)
phi_v = m.Const(virus.phi)
xi_v = m.Const(virus.xi)

SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))

m.Equation(x1_h==x0_h.dt())
m.Equation(x2_h==x1_h.dt())
m.Equation(x2_h == ((SU0 / (ST0**2))*alpha_h*u2 + (SU0 / ST0)*beta_h*u1 + gamma_h*SU0*u0 - phi_h*(SX0 / ST0)*x1_h - xi_h*SX0*x0_h)/(SX0 / (ST0**2)))

m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True) # solve



first_derivative = (np.diff(control_input)/dt)
second_derivative = np.diff(first_derivative)/dt

t_ = t[:-2]
control_input = control_input[:-2]
first_derivative = first_derivative[:-1]

sum = host_cell.alpha * ((T0*T0))*second_derivative + host_cell.beta * ((T0))*first_derivative + host_cell.gamma * control_input


u1 = np.loadtxt( 'PSOPT/build/u1.dat', dtype=np.float128)[:-1]
u2 = np.loadtxt( 'PSOPT/build/u2.dat', dtype=np.float128 )[:-1]
sim_x0_v = np.loadtxt( 'PSOPT/build/x0_v.dat' )[:-1]
sim_x0_h = np.loadtxt( 'PSOPT/build/x0_h.dat' )[:-1]


plt.figure(1)
plt.plot(t_,host_cell.gamma * control_input,marker='o')
plt.plot(sim_t,host_cell.gamma * sim_values)
plt.show()
plt.figure(2)
plt.plot(sim_t, u1, marker='o')
plt.plot(t_, first_derivative*T0,marker='o')
plt.show()
plt.figure(3)
plt.plot(sim_t, u2, marker='o')
plt.plot(t_, second_derivative*T0*T0,marker='o')
plt.show()
plt.figure(4)
plt.plot(t_, sum,marker='o')
plt.show()
plt.figure(5)
plt.plot(t, convolve_output(control_input, host_cell, dt))
plt.plot(t, convolve_output(control_input, virus, dt))
plt.plot(sim_t, sim_x0_v/1e6)
plt.plot(sim_t, sim_x0_h/1e6)
plt.figure(6)
plt.plot(m.time*T0, x0_v.value)
plt.plot(m.time*T0, x0_h.value)
plt.show()



print((np.max(convolve_output(ideal_values, host_cell, dt) * 1e6) - np.min(convolve_output(ideal_values, host_cell, dt) * 1e6))
            /(np.max(convolve_output(ideal_values, virus, dt) * 1e6) - np.min(convolve_output(ideal_values, virus, dt) * 1e6)))






#
