from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO
import numpy as np
from scipy import signal
from scipy.ndimage.interpolation import shift
from scipy import interpolate


m = GEKKO() # initialize gekko
m.options.MAX_ITER = 30
nt = 1000
end = 1e-8

T0 = 1e-8
U0 = 1.0
X0 = 1.0
m.time = np.linspace(0, end/T0, nt, dtype=np.float128)
t1 = np.linspace(0, end, nt, dtype=np.float128)



x0_v = m.Var(value=0)
x1_v = m.Var(value=0)
x2_v = m.Var(value=0)

t = m.Param(value=m.time)

sim_t = (np.loadtxt( 'PSOPT/build/t.dat', dtype=np.float128) * T0)[:-1]
sim_u0 = np.loadtxt( 'PSOPT/build/u0.dat', dtype=np.float128 )[:-1]

interpolant = interpolate.InterpolatedUnivariateSpline(sim_t, sim_u0, k=5)

# host_cell = 
virus = default_host_cell(t1)



u_step = interpolant(np.float64(m.time)*T0)
plt.plot(t1, u_step)
plt.show()
u0 = m.Param(value=u_step)
u1_ = np.diff(u_step, prepend=u_step[0])/((end/nt)/T0)
u1 = m.Param(value=u1_)
u2_ = np.diff(u1_, prepend=u1_[0])/((end/nt)/T0)
u2 = m.Param(value=u2_)

alpha_v = m.Const(virus.alpha)
beta_v = m.Const(virus.beta)
gamma_v = m.Const(virus.gamma)
phi_v = m.Const(virus.phi)
xi_v = m.Const(virus.xi)

SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))

m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True)

virus_output = np.array(x0_v.value)

plt.figure(1) # plot results
integration_method = virus_output*X0 / U0

plt.plot(m.time*T0,integration_method,'r',marker='o', label="GEKKO integrator")

convolution_method_output = convolve_output(u_step, virus, end/nt)
plt.plot(t1, convolution_method_output, label="lib_onvolution_method_output")



system = signal.lti(np.float64([virus.a_3,virus.a_2,virus.a_1]),np.float64([virus.b_3,virus.b_2,virus.b_1]))
#output is transfer function times input times R!

#odeint did not produce the correct result except when rtol and atol were
_, yout, _ = signal.lsim2(system, U=np.float64(u_step), T=np.float64(t1),rtol=1e-11,atol=1e-11)
plt.plot(t1, yout*virus.R, label='lsim2')
_, yout, _ = signal.lsim(system, U=np.float64(u_step), T=np.float64(t1))
plt.plot(t1, yout*virus.R, label='lsim')

plt.legend()

plt.show()
