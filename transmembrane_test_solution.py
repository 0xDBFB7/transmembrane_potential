from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os
from scipy import interpolate

T0 = 1e-8


sim_t = np.loadtxt( 'PSOPT/build/t.dat' ) * T0


t0 = 0
tstop = sim_t[-1]
dt = 0.0005e-9
t = np.linspace(t0, tstop, int(tstop/dt))

# dt = (t[1] - t[0])

sim_values = np.loadtxt( 'PSOPT/build/u0.dat' )


# flinear = interpolate.interp1d(x, sim_values)
fcubic = interpolate.interp1d(sim_t, sim_values, kind='linear')

ideal_values = fcubic(t)

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)

# small_host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 2.5e-6, 5e-9, t)
# host_cell = small_host_cell

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)

def normalized_gaussian_pulse(t,fwhm):
    sigma = fwhm/2.355
    return np.exp(-((t**2.0)/(2.0*(sigma**2.0))))


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
