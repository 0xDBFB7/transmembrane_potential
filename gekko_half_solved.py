from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt
from transmembrane_lib import *

t = np.array([0])

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 50e-6, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)


'''
Kwin van der Veen on math.se clued me in:

Without additional constraints or penalty to u(t) I think that this problem is not well posed.
Namely, one could use a solution such that Ahu¨(t)+Bhu˙(t)+Chu(t)=0 and scale it by any arbitrary large scalar.


that is to say, every solution to xh==
'''



m = GEKKO() # initialize gekko
m.options.MAX_ITER = 200
nt = 800

T0 = 1e-7
U0 = 1.0
X0 = 1e-6
end = 1e-6 / T0

m.time = np.linspace(0,end,nt)

# Variables
x0_v = m.Var(value=0)
x1_v = m.Var(value=0)
x2_v = m.Var(value=0)

t = m.Param(value=m.time)

u0 = m.Var(value=0)

# m.Equation(u0 == 1) #doesn't seem to behave well with this discontinuity.
# the gaussian pulse below works much better.
# m.Equation(u0 == m.sin(t)) # for simulation
# m.Equation(u0 == m.exp(-((((t*T0)-((end*T0))/2.0))**2.0)/(2.0*(((((end*T0)/10.0))**2.0))))) # for simulation

u1 = m.Var(value=1)
m.Equation(u1==u0.dt())
u2 = m.Var()
m.Equation(u2==u1.dt())

# goes here!

SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))


m.Equation(0 == ((SU0 / (ST0**2))*alpha_h*u2 + (SU0 / ST0)*beta_h*u1 + gamma_h*SU0*u0)/(SX0 / (ST0**2)))

#0 = (𝛂ₕ*d²u/dt² + 𝛃ₕ*du/dt + 𝛄ₕ*u)


m.Obj(-m.integral(x0_v*x0_v))
# m.options.OTOL = 1e-6
# m.options.RTOL = 1e-6

m.options.IMODE = 6 # optimal control mode

# m.options.IMODE = 4 # dynamic simulation

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
plt.figure(2) # plot results
plt.plot(m.time,u0.value,'g',label=r'$u$')


plt.show()
