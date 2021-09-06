from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os

FWHM = 0.5e-7
t0 = 0
tstop = FWHM*8
dt = 0.005e-8
t = np.linspace(t0, tstop, int(tstop/dt))


def normalized_gaussian_pulse(t,fwhm):
    sigma = fwhm/2.355
    return np.exp(-((t**2.0)/(2.0*(sigma**2.0))))

input = np.zeros_like(t)

ideal_values = normalized_gaussian_pulse(t-tstop/2.0,FWHM)


host_cell = Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)
virus = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8), np.float128(60), np.float128(50e-9), np.float128(14e-9),t)

host_cell_permeabilized = Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7 + 0.2), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)
virus_permeabilized = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8 + 0.2), np.float128(60), np.float128(50e-9), np.float128(14e-9),t)

plt.style.use('grayscale')
# plt.plot(t/1e-9, ideal_values, linestyle="solid", label="Input values")
host_output = convolve_output(ideal_values, host_cell, dt) * 1e6
virus_output = convolve_output(ideal_values, virus, dt) * 1e7
host_ep_output = convolve_output(ideal_values, host_cell_permeabilized, dt) * 1e7
virus_ep_output = convolve_output(ideal_values, virus_permeabilized, dt) * 1e7

plt.plot(t/1e-9, host_output, linestyle="solid", label="Intact host cell scaled 1/10 ($\sigma_{membrane} = 10^-7$)")
plt.plot(t/1e-9, virus_output, linestyle="dashed", label="Intact virus ($\sigma_{membrane} = 10^-8$)")
plt.plot(t/1e-9, host_ep_output, linestyle="dashdot", label="Host cell ($\sigma_{membrane} = 10^-7$)")
plt.plot(t/1e-9, virus_ep_output, linestyle="dotted", label="Virus perm")

plt.xlabel("Time (nanoseconds)\nY axis: 1 V/m input waveform, membrane $\Delta$V scaled by $10^7$")
plt.legend()
plt.savefig("plots/pre_coefficient_switchover.svg")

plt.clf()

host_cell = Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)
virus = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8), np.float128(5), np.float128(50e-9), np.float128(5e-9),t)

host_cell_permeabilized = Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7 + 0.2), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)
virus_permeabilized = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8 + 0.2), np.float128(5), np.float128(50e-9), np.float128(5e-9),t)

plt.style.use('grayscale')
# plt.plot(t/1e-9, ideal_values, linestyle="solid", label="Input values")
host_output = convolve_output(ideal_values, host_cell, dt) * 1e6
virus_output = convolve_output(ideal_values, virus, dt) * 1e7
host_ep_output = convolve_output(ideal_values, host_cell_permeabilized, dt) * 1e7
virus_ep_output = convolve_output(ideal_values, virus_permeabilized, dt) * 1e7

plt.plot(t/1e-9, host_output, linestyle="solid", label="Intact host cell scaled 1/10 ($\sigma_{membrane} = 10^-7$)")
plt.plot(t/1e-9, virus_output, linestyle="dashed", label="Intact virus ($\sigma_{membrane} = 10^-8$)")
plt.plot(t/1e-9, host_ep_output, linestyle="dashdot", label="Host cell ($\sigma_{membrane} = 10^-7$)")
plt.plot(t/1e-9, virus_ep_output, linestyle="dotted", label="Virus perm")

plt.xlabel("Time (nanoseconds)\nY axis: 1 V/m input waveform, membrane $\Delta$V scaled by $10^7$")
plt.legend()
plt.savefig("plots/post_coefficient_switchover.svg")


plt.show()




#
