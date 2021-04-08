from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle

t0 = 0
tstop = 50e-6
dt = 100e-9
t = np.linspace(t0, tstop, int(tstop/dt))

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)

input = np.ones_like(t)
input[0:(len(input)//100)] = 0
input[(len(input)//2):-1] = 0


output = convolve_output(input, host_cell, dt) * 1e6

# plt.plot(np.arange(len(output)) * dt, output)
# print(total_waveform_energy(output,dt))

# plt.plot(t, input)


def cost_function(input):
    print(np.max(input))
    host_cell_output = total_waveform_energy(convolve_output(input, host_cell, dt) * 1e6, dt)
    virus_output = total_waveform_energy(convolve_output(input, virus, dt) * 1e6, dt)
    print(host_cell_output, virus_output)
    # plt.plot(np.arange(len(output)) * dt,input)
    # plt.plot(np.arange(len(output)) * dt,convolve_output(input, host_cell, dt) * 1e6)
    # plt.show()
    return (host_cell_output / virus_output )


# x0 = np.ones_like(t)
x0 = input

tubthumper = basinhopping
minimizer_kwargs = dict(method="L-BFGS-B", options={"disp":True})



filename = 'globalsave.pkl'
try:
    dill.load_session(filename)
except:
    ideal_values = tubthumper(cost_function, x0, niter=100, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=5)["x"]
    dill.dump_session(filename)


# you may not like it, but this is the ideal

plt.plot(t, ideal_values)
plt.plot(t, convolve_output(ideal_values, host_cell, dt) * 1e6)
plt.plot(t, convolve_output(ideal_values, virus, dt) * 1e6)
plt.show()



# x = minimize(cost_function, x0, options={"disp":True}).x
