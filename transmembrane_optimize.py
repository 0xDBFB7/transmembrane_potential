from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os

t0 = 0
tstop = 1e-6
dt = 2e-9
t = np.linspace(t0, tstop, int(tstop/dt))

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 80, 50e-9, 14e-9, t)

print(f"Host cell {host_cell.tau_1}, {host_cell.tau_2}")
print(f"virus {virus.tau_1}, {virus.tau_2}")

#20e-6/50e-9

#input = np.ones_like(t)
input = np.sin(t)
input[0:(len(input)//100)] = 0
input[(len(input)//2):-1] = 0
input[-1] = 0


output = convolve_output(input, host_cell, dt) * 1e6

# plt.plot(np.arange(len(output)) * dt, output)
# print(total_waveform_energy(output,dt))

# plt.plot(t, input)


def cost_function(input):
    # print(np.max(input))

    host_cell_output =np.linalg.norm(convolve_output(input, host_cell, dt) * 1e6)
    # host_cell_output = np.sum(host_cell_output[host_cell_output > 0])
    virus_output = np.linalg.norm(convolve_output(input, virus, dt) * 1e6)
    # virus_output = np.sum(virus_output[virus_output > 0])


    # print(host_cell_output, virus_output)
    # plt.plot(np.arange(len(output)) * dt,input)
    # plt.plot(np.arange(len(output)) * dt,convolve_output(input, host_cell, dt) * 1e6)
    # plt.show()

    #tstop -
    return ((host_cell_output - virus_output) + (abs(0.0 - np.linalg.norm(input)))) # put a total energy limit here?

def diagnostics(xk):
    input = xk
    os.system("clear")
    host_cell_output =np.linalg.norm(convolve_output(input, host_cell, dt) * 1e6)
    # host_cell_output = np.sum(host_cell_output[host_cell_output > 0])
    virus_output = np.linalg.norm(convolve_output(input, virus, dt) * 1e6)
    print((host_cell_output / virus_output))
    print(np.linalg.norm(input))
# x0 = np.ones_like(t)
x0 = input

tubthumper = basinhopping
minimizer_kwargs = dict(method="L-BFGS-B", options={"disp":True})

#wait, this is a inear combination except for the two time constants, right?
#so how could there ever be a difference?
#
filename = 'globalsave.pkl'
try:
    dill.load_session(filename)
    print("LOADED PREVIOUS SESSION")
except:
    ideal_values = minimize(cost_function, x0, method="BFGS", options={"disp":True}, callback=diagnostics).x
    #ideal_values = tubthumper(cost_function, x0, niter=100, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=5)["x"]
    dill.dump_session(filename)

# you may not like it, but this is the ideal

plt.plot(t, ideal_values)
plt.plot(t, convolve_output(ideal_values, host_cell, dt) * 1e6)
plt.plot(t, convolve_output(ideal_values, virus, dt) * 1e6)

plt.show()




#
