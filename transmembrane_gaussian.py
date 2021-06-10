from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os

FWHM = 0.5e-9
t0 = 0
tstop = FWHM*8
dt = 0.005e-9
t = np.linspace(t0, tstop, int(tstop/dt))

host_cell = default_host_cell(t)
virus = default_virus(t)

print(epsilon_0)

print(f"Host cell {host_cell.tau_1}, {host_cell.tau_2}")
print(f"virus {virus.tau_1}, {virus.tau_2}")

#20e-6/50e-9

def normalized_gaussian_pulse(t,fwhm):
    sigma = fwhm/2.355
    return np.exp(-((t**2.0)/(2.0*(sigma**2.0))))




input = np.zeros_like(t)

ideal_values = normalized_gaussian_pulse(t-tstop/2.0 - FWHM,FWHM) - normalized_gaussian_pulse(t-tstop/2.0 + FWHM,FWHM)

plt.style.use('grayscale')
plt.plot(t, ideal_values, linestyle="solid", label="Input values")
host_output = convolve_output(ideal_values, host_cell, dt) * 1e7
virus_output = convolve_output(ideal_values, virus, dt) * 1e7
plt.plot(t, host_output, linestyle="dotted", label="Host cell")
plt.plot(t, virus_output, linestyle="dashed", label="Virus")

plt.xlabel("Time (nanoseconds)")
plt.ylabel("1 V/m input waveform,\nTransmembrane voltage scaled by $10^7$")
plt.legend()

plt.savefig("plots/gaussian_pulse.png")



print(np.max(np.abs(host_output))/np.max(np.abs(virus_output)))
print(np.sum(host_output**2.0) / np.sum(virus_output**2.0))

plt.show()




#
