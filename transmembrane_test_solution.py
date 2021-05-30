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
sim_u0 = np.loadtxt( 'PSOPT/build/u0.dat', dtype=np.float128 )[:-1]

# t0 = 0
# tstop = sim_t[-1]
# print(tstop)
# dt = 0.0005e-9 / T0
# t = np.linspace(t0, tstop, int(tstop/dt), dtype=np.float128)

interpolant = interpolate.InterpolatedUnivariateSpline(sim_t, sim_values, k=5)
control_input = interpolant(np.float64(t))

host_cell = default_host_cell(sim_t * T0)
virus = default_virus(sim_t * T0)





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

#
# plt.figure(1)
# plt.plot(t_,host_cell.gamma * control_input,marker='o')
# plt.plot(sim_t,host_cell.gamma * sim_values)
# plt.show()
# plt.figure(2)
# plt.plot(sim_t, u1, marker='o')
# plt.plot(t_, first_derivative*T0,marker='o')
# plt.show()
# plt.figure(3)
# plt.plot(sim_t, u2, marker='o')
# plt.plot(t_, second_derivative*T0*T0,marker='o')
# plt.show()
# plt.figure(4)
# plt.plot(t_, sum,marker='o')
# plt.show()
# plt.figure(5)
# plt.plot(t, convolve_output(control_input, host_cell, dt))
# plt.plot(t, convolve_output(control_input, virus, dt))
# plt.plot(sim_t, sim_x0_v/1e6)
# plt.plot(sim_t, sim_x0_h/1e6)
plt.figure(6)
plt.plot(m.time*T0, x0_v.value)
plt.plot(m.time*T0, x0_h.value)
plt.show()






#
