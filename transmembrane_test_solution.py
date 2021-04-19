from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os

# t0 = 0
# tstop = 1e-6
# dt = 0.005e-9
# t = np.linspace(t0, tstop, int(tstop/dt))

t = np.loadtxt( 'PSOPT/build/t.dat' )
dt = (t[1] - t[0])

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)

# small_host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 2.5e-6, 5e-9, t)
# host_cell = small_host_cell

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)

def normalized_gaussian_pulse(t,fwhm):
    sigma = fwhm/2.355
    return np.exp(-((t**2.0)/(2.0*(sigma**2.0))))

ideal_values = np.loadtxt( 'PSOPT/build/u0.dat' )

plt.figure(1)
plt.plot(t, ideal_values)
plt.show()
plt.figure(2)
plt.plot(t, convolve_output(ideal_values, host_cell, dt))
plt.plot(t, convolve_output(ideal_values, virus, dt))

print((np.max(convolve_output(ideal_values, host_cell, dt) * 1e6) - np.min(convolve_output(ideal_values, host_cell, dt) * 1e6))
            /(np.max(convolve_output(ideal_values, virus, dt) * 1e6) - np.min(convolve_output(ideal_values, virus, dt) * 1e6)))


plt.show()




#
