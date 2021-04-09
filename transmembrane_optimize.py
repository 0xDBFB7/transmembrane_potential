from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle

t0 = 0
tstop = 2e-6
dt = 0.005e-9
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


# def cost_function(input):
#     # print(np.max(input))
#
#     host_cell_output = np.sum(convolve_output(input, host_cell, dt) * 1e6)
#     virus_output = np.sum(convolve_output(input, virus, dt) * 1e6)
#
#
#     # print(host_cell_output, virus_output)
#     # plt.plot(np.arange(len(output)) * dt,input)
#     # plt.plot(np.arange(len(output)) * dt,convolve_output(input, host_cell, dt) * 1e6)
#     # plt.show()
#
#
#     print((host_cell_output / virus_output))
#     return (host_cell_output / virus_output) + (1.0 - total_waveform_energy(input,dt)) # put a total energy limit here?

# def cost_function(input):

    # return (host_cell_output / virus_output)) # put a total energy limit here?

# x0 = np.ones_like(t)
x0 = input

tubthumper = basinhopping
minimizer_kwargs = dict(method="L-BFGS-B", options={"disp":True})

# take the FT (fft?), put as much power into that spectrum as possible?

#wait, this is a inear combination except for the two time constants, right?
#so how could there ever be a difference?
#
# filename = 'globalsave.pkl'
# try:
#     dill.load_session(filename)
#     print("LOADED PREVIOUS SESSION")
# except:
#     ideal_values = minimize(cost_function, x0, method="CG", options={"disp":True}).x
#     #ideal_values = tubthumper(cost_function, x0, niter=100, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=5)["x"]
#     dill.dump_session(filename)

# you may not like it, but this is the ideal

# plt.plot(t, ideal_values)
# plt.plot(t, convolve_output(ideal_values, host_cell, dt) * 1e6)
# plt.plot(t, convolve_output(ideal_values, virus, dt) * 1e6)
#
# plt.show()


plt.figure()
plt.plot(t, host_cell.step_response / np.linalg.norm(host_cell.step_response))
plt.plot(t, virus.step_response / np.linalg.norm(virus.step_response))
plt.show()

n = 1000000
host_cell_fft = np.fft.fft(host_cell.step_response / np.linalg.norm(host_cell.step_response),n)
virus_fft = np.fft.fft(virus.step_response / np.linalg.norm(virus.step_response),n)
freq = np.fft.fftfreq(n,dt)

plt.figure()
plt.plot(freq, virus_fft-host_cell_fft)
plt.show()
plt.figure()

shaped_pulse = np.fft.ifft(virus_fft-host_cell_fft)
shaped_dt = dt*(len(t) / len(shaped_pulse)) # to get high freq res, need a large fft n; but this messes up the timestep
# shaped_pulse = shaped_pulse[0:len(shaped_pulse) //3]
shaped_pulse -= np.min(np.real(shaped_pulse)) #add a DC offset
shaped_pulse /= np.max(shaped_pulse) #normalize
# shaped_pulse = np.pad(shaped_pulse,[0,500000])
# shaped_pulse = np.tile(shaped_pulse, 20)

shaped_times = np.arange(len(shaped_pulse)) * shaped_dt

plt.plot(shaped_times,np.real(shaped_pulse))
plt.show()

plt.figure()
plt.plot(shaped_times,convolve_output(shaped_pulse, host_cell, shaped_dt))
plt.plot(shaped_times,convolve_output(shaped_pulse, virus, shaped_dt))
plt.show()



#
