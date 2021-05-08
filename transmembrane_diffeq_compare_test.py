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


host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, t1 )
virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t1 )

 # differential equation is wrong - must be just after gekko_.

control_input = np.sin(t1/((end/30.0)/T0))



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

alpha_h = m.Const(host_cell.R*host_cell.a_1/host_cell.b_1)
beta_h = m.Const(host_cell.R*host_cell.a_2/host_cell.b_1)
gamma_h = m.Const(host_cell.R*host_cell.a_3/host_cell.b_1)
phi_h = m.Const(host_cell.b_2/host_cell.b_1)
xi_h = m.Const(host_cell.b_3/host_cell.b_1)

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
plt.plot(m.time*T0,virus_output*X0 / U0,'r',label=r'$x0_v$')
# plt.plot(m.time*T0,host_output*X0 / U0,'b',label=r'$x0_h$')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()
plt.figure(2) # plot results
plt.plot( t1, convolve_output(control_input, virus, end/nt))

print(np.max((virus_output*X0 / U0) - control_input))


plt.show()
