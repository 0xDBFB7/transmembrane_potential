from transmembrane_lib import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize,basinhopping
import dill
from pytexit import py2tex
import pickle
import os

t0 = 0
tstop = 2e-6
dt = 0.01e-9
t = np.linspace(t0, tstop, int(tstop/dt))

# host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t)
#
# virus = Cell(0.3, 80, 0.005, 30, 1e-8, 80, 50e-9, 14e-9, t)

host_cell = default_host_cell(t)

virus = default_virus(t)

print(f"Host cell {host_cell.tau_1 / 1e-9}, {host_cell.tau_2 / 1e-9}")
print(f"virus {virus.tau_1 / 1e-9}, {virus.tau_2 / 1e-9}")

ideal_values = np.ones_like(t)
ideal_values[0:20] = 0

host_output = convolve_output(ideal_values, host_cell, dt)
virus_output = convolve_output(ideal_values, virus, dt)

plt.style.use('grayscale')
plt.plot(t/1e-9, host_output / np.max(host_output),linestyle="dashed", label="Host cell")
plt.plot(t/1e-9, virus_output / np.max(virus_output), label="Virus")

plt.xlabel("Time (nanoseconds)")
plt.ylabel("Normalized transmembrane voltage")
plt.legend()

plt.savefig("plots/step_course.png")

plt.show()




#
