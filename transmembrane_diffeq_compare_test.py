from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from gekko import GEKKO
import numpy as np



m = GEKKO() # initialize gekko
m.options.MAX_ITER = 100000
nt = 301
end = 1e-6


t0 = 0
tstop = 1e-6

T0 = 1e-6
U0 = 1.0
X0 = 1e-6
m.time = np.linspace(t0, tstop/T0, nt)
t1 = np.linspace(t0, tstop, nt)


virus = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8), np.float128(60), np.float128(50e-9), np.float128(14e-9), t1)

 # differential equation is wrong - must be just after gekko_.

control_input = np.sin(t1/((end/30.0)))



# Variables
x0_v = m.Var(value=0)
x1_v = m.Var()
x2_v = m.Var()

t = m.Param(value=m.time)

u0 = m.Var()
# m.Equation(u0 == 1)
# the gaussian pulse below works much better.
m.Equation(u0 == m.sin(t/((end/30)/T0))) # for simulation
# m.Equation(u0 == m.exp(-((((t*T0)-((end*T0))/2.0))**2.0)/(2.0*(((((end*T0)/10.0))**2.0))))) # for simulation

u1 = m.Var()
m.Equation(u1==u0.dt())
u2 = m.Var()
m.Equation(u2==u1.dt())

alpha_v = m.Const(virus.R*virus.a_1/virus.b_1)
beta_v = m.Const(virus.R*virus.a_2/virus.b_1)
gamma_v = m.Const(virus.R*virus.a_3/virus.b_1)
phi_v = m.Const(virus.b_2/virus.b_1)
xi_v = m.Const(virus.b_3/virus.b_1)


SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))



m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True) # solve


virus_output = np.array(x0_v.value)
# host_output = np.array(x0_h.value)

# print(np.max(virus_output) / np.max(host_output))

plt.figure(1) # plot results
# plt.plot(m.time,x1.value,'k-',label=r'$x_1$')
# plt.plot(m.time,x2.value,'b-',label=r'$x_2$')
integration_method = virus_output*X0 / U0
plt.plot(m.time*T0,integration_method,'r',label=r'$x0_v$')
# convolution_method_output = convolve_output(control_input, virus, end/nt)
# plt.plot( t1, convolution_method_output)


print(np.mean(np.abs(integration_method - convolution_method_output))/(np.max(convolution_method_output)))


plt.show()
