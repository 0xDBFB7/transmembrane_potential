from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np

t0 = 0
tstop = 10.0
t = np.linspace(t0, tstop, 1000)

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9)

# deltaV = delta_transmembrane_unit_ramp(t, host_cell)

#asymmetric rise/fall trapezoid
deltaV = delta_transmembrane_trapezoid(t, 1, 2.0, 2.0, 1.0, host_cell, 1.0)

#square wave
# deltaV = delta_transmembrane_trapezoid(t, 1, 1e-9, 1e-9, 1.0, host_cell, 1.0)



plt.plot(t,deltaV)
plt.show()
