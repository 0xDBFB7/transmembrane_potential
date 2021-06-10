from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO
import numpy as np
from scipy import signal
from scipy.ndimage.interpolation import shift
from scipy import interpolate
from scipy.optimize import differential_evolution, Bounds, NonlinearConstraint
end = 1e-7
nt = 100
t = np.linspace(0, end, nt, dtype=np.float64)#float128

virus = default_virus(np.array([]))
host_cell = default_host_cell(np.array([]))


virus_system = signal.lti(np.float64([virus.a_3,virus.a_2,virus.a_1]),np.float64([virus.b_3,virus.b_2,virus.b_1]))
host_system = signal.lti(np.float64([host_cell.a_3,host_cell.a_2,host_cell.a_1]), np.float64([host_cell.b_3,host_cell.b_2,host_cell.b_1]))

def cost_function(control_input):
    #lsim2 uses ODE integration,
    #lsim uses convolution (maybe FFT in high N cases, haven't checked)
    _, x_v, _ = signal.lsim2(virus_system, U=control_input, T=t, atol=1e-9)
    _, x_h, _ = signal.lsim2(host_system, U=control_input, T=t, atol=1e-9)#rtol e-11
    rho = 1.0
    return -rho*np.sum(x_v*x_v) + np.sum(x_h*x_h)

def constr_f(x):
    return np.sum(x)

nlc = NonlinearConstraint(constr_f, 0.001, 0.1)

# def callback(x,convergence=None):
#     plt.clf()
#     plt.plot(t,x_v)
#     plt.plot(t,x_h)
#     plt.draw()
#     plt.pause(0.001)

lb = np.ones_like(t)*-1.0
lb[0] = 0.0
ub = np.ones_like(t)
ub[0] = 0.0
bound = Bounds(lb,ub)

result = differential_evolution(cost_function,bound,maxiter=2,constraints=(nlc),workers=14,disp=True)["x"] #callback=callback,

plt.figure(0)
plt.plot(t,result)
_, x_v, _ = signal.lsim2(virus_system, U=result, T=t)
_, x_h, _ = signal.lsim2(host_system, U=result, T=t)#, rtol=1e-8, atol=1e-8
plt.figure(1)
plt.plot(t,x_v)
plt.plot(t,x_h)
plt.legend()
plt.show()
