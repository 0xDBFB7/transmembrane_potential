from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO
import numpy as np
from scipy import signal
from scipy.ndimage.interpolation import shift
from scipy import interpolate

end = 1e-8
nt = 3000
t = np.linspace(0, end, nt, dtype=np.float128)

virus = default_virus(t)
host_cell = default_host_cell(t)


system = signal.lti(np.float64([virus.a_3,virus.a_2,virus.a_1]),np.float64([virus.b_3,virus.b_2,virus.b_1]))

_, yout, _ = signal.lsim2(system, U=np.float64(u_step), T=np.float64(t),rtol=1e-11,atol=1e-11)


plt.legend()
plt.show()
