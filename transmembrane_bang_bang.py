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
dt = 0.1e-9
t = np.linspace(t0, tstop, int(tstop/dt))

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 80, 50e-9, 14e-9, t)

print(f"Host cell {host_cell.tau_1}, {host_cell.tau_2}")
print(f"virus {virus.tau_1}, {virus.tau_2}")

#20e-6/50e-9

input = np.zeros_like(t)
# input = np.sin(t)
# input[0:(len(input)//100)] = 0
# input[(len(input)//2):-1] = 0
# input[-1] = 0


# output = convolve_output(input, host_cell, dt) * 1e6

# plt.plot(np.arange(len(output)) * dt, output)
# print(total_waveform_energy(output,dt))

# plt.plot(t, input)
output = np.zeros_like(t)



for i in range(0, len(t)):

    options = np.array([0.00001, -1, 1])
    ratios = np.zeros(3)

    for op_idx,opt in enumerate(options):
        input[i] = opt
        host_cell_output = abs((convolve_output(input, host_cell, dt) * 1e6)[i])
        virus_output = abs((convolve_output(input, virus, dt) * 1e6)[i])
        # if(host_cell_output

        # ratios[op_idx] = virus_output / host_cell_output

    print(ratios)

    input[i] = options[np.argsort(ratios)][-1]


# plt.plot(t, input)

    # output = [i]

    # print(host_cell_output, virus_output)
    # plt.plot(np.arange(len(output)) * dt,input)
    # plt.plot(np.arange(len(output)) * dt,convolve_output(input, host_cell, dt) * 1e6)
    # plt.show()

    #tstop -
    # return host_cell_output - virus_output
    # return ((host_cell_output - virus_output) + (abs(1.0 - np.linalg.norm(input))))
    # return ((host_cell_output - virus_output) + (abs(1.0 - np.max(np.abs(input)))))

#I think something's wrong with the convolve?
# not obeying causality?




# x0 = np.ones_like(t)


#wait, this is a inear combination except for the two time constants, right?
#so how could there ever be a difference?
#
# filename = 'globalsave.pkl'
# try:
#     dill.load_session(filename)
#     print("LOADED PREVIOUS SESSION")
# except:
#     ideal_values = minimize(cost_function, x0, method="BFGS", options={"disp":True}, callback=diagnostics).x
#     #gtol:
#     #ideal_values = tubthumper(cost_function, x0, niter=100, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=5)["x"]
#     dill.dump_session(filename)

# you may not like it, but this is the ideal

ideal_values = input

plt.plot(t, ideal_values)
plt.plot(t, convolve_output(ideal_values, host_cell, dt) * 1e6)
plt.plot(t, convolve_output(ideal_values, virus, dt) * 1e6)



plt.show()




#
