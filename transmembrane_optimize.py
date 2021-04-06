from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping

t0 = 0
tstop = 100e-6
t = np.linspace(t0, tstop, 100000)

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9)


# deltaV = delta_transmembrane_unit_ramp(t, host_cell)

#asymmetric rise/fall trapezoid
# deltaV = delta_transmembrane_trapezoid(t, 1, 2.0, 3.0, 1.0, host_cell, 1.0)

#square wave
# deltaV = delta_transmembrane_trapezoid(t, 100e-6, 2.3e-9, 2.3e-9, 30e-6, host_cell, 10.0e6)
# plt.plot(t,deltaV)
# deltaV = delta_transmembrane_rectangular_train(t, 1e-6, 3e-6, 2, 10e-6, host_cell, 10.0e6)
# plt.plot(t, deltaV)

def cost_function(x):
    # x = np.abs(x)
    # print(x)
    host_amp = np.max(delta_transmembrane_rectangular_array(t, 10, x, host_cell))
    vir_amp = np.max(delta_transmembrane_rectangular_array(t, 10, x, virus))
    return host_amp - vir_amp
    # return np.max(delta_transmembrane_rectangular_train(t, 1e-6, x[0], 10, x[1], host_cell, 10.0e6)) - np.max(delta_transmembrane_rectangular_train(t, 1e-6, x[0], 10, x[1], virus, 10.0e6))

#use rectangle pulses to construct a waveform?
#
# x0 = [1e-9, 1e-6]
x0 = np.zeros(30)
x = minimize(cost_function, x0, method='CG', options={"disp":True}).x
#x= [1.49011602e-08, 0.1e-06]
#still need polarities
#x=[1.17462645e-04,  1.03637294e-06]
deltaV = delta_transmembrane_rectangular_train(t, 1e-6, x[0], 5, x[1], host_cell, 10.0e6)
plt.plot(t, deltaV)
deltaV = delta_transmembrane_rectangular_train(t, 1e-6, x[0], 5, x[1], virus, 10.0e6)
plt.plot(t, deltaV)

# deltaV =   delta_transmembrane_trapezoid(t, 0e-6, 2.3e-9, 2.3e-9, 300000e-6, virus, 100.0e6)
# deltaV +=  delta_transmembrane_trapezoid(t, 0e-6, 2.3e-9, 2.3e-9, 300000e-6, virus, 100.0e6)
# plt.plot(t, deltaV)


# plt.plot(t,deltaV)


plt.show()
