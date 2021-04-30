from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

from transmembrane_lib import *

t = np.array([0])

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 50e-6, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)



m = GEKKO() # initialize gekko
m.options.MAX_ITER = 100000
nt = 401

X0 = 1.0 / (virus.b_3/virus.b_1)
U0 = 1.0 / (virus.R*virus.a_3/virus.b_1)
T0 = 1.0 / (virus.a_3 / virus.a_2)



end = 1e-6 / T0
m.time = np.linspace(0,end,nt)

# Variables
x0_v = m.Var(value=0)
x1_v = m.Var(value=0)
x2_v = m.Var(value=0)

# x0_h = m.Var(value=0)
# x1_h = m.Var(value=0)
# x2_h = m.Var(value=0)

t = m.Param(value=m.time)



SX0 = m.Const(1.0 / (virus.b_3/virus.b_1))
ST0 = m.Const(1.0 / (virus.a_3 / virus.a_2))

print(end)

u0 = m.Var(value=1)

m.Equation(u0 == 1)
# m.Equation(u0 == m.sin(t)) # for simulation
# m.Equation(u0 == m.exp(-(((t-(end/2.0))**2.0)/(2.0*(((0.1e-7)/t_0)**2.0))))) # for simulation


u1 = m.Var()
m.Equation(u1==u0.dt())
u2 = m.Var()
m.Equation(u2==u1.dt())


a1_v = m.Const(virus.a_1)
a2_v = m.Const(virus.a_2)
a3_v = m.Const(virus.a_3)
b1_v = m.Const(virus.b_1)
b2_v = m.Const(virus.b_2)
b3_v = m.Const(virus.b_3)
R_v = m.Const(virus.R)


# a1_h = m.Const(host_cell.a_1)
# a2_h = m.Const(host_cell.a_2)
# a3_h = m.Const(host_cell.a_3)
# b1_h = m.Const(host_cell.b_1)
# b2_h = m.Const(host_cell.b_2)
# b3_h = m.Const(host_cell.b_3)
# R_h = m.Const(host_cell.R)


# alpha = ((a1_h*a3_v*b1_v*R_h)/((a2_v*a2_v)*b1_h*R_v))
# beta = (((1/(R_v*(a3_v/b1_v)))/(1 / (a3_v / a2_v)))*(R_h*(a2_h/b1_h)))
# gamma = ((a3_h*b1_v*R_h)/(a3_v*b1_h*R_v))
# psi = ((a2_v*b2_h)/(a3_v*b1_h))
# xi = (((a2_v*a2_v)*b3_h)/((a3_v*a3_v)*b1_h))


print((((a2_v*a2_v)*b3_v)/((a3_v*a3_v)*b1_v)))

# Equations

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == (((a1_v*a3_v)/(a2_v*a2_v)) * u2 + u1 + u0 - ((a3_v*b2_v)/(a2_v*b3_v))*x1_v - x0_v) / (SX0 / ((ST0)**2.0)))





# m.Equation(x1_h==x0_h.dt())
# m.Equation(x2_h==x1_h.dt())
# m.Equation(x2_h == alpha*u2 + beta*u1 + gamma*u0 - psi*x1_h - xi*x0_h)

#
# int_h = m.Var()
# m.Equation(int_h==m.integral(x0_h**2.0))
# m.Equation(t==m.vsum(m.abs2(u0)))

# m.Equation(m.integral(u0*u0)==(end))
# m.Equation(m.vsum(m.abs2(u0))==0.1)
# m.Equation(m.integral(m.abs2(u0))==0.1)
# integral()
# abs2()


# m.Obj(m.integral(x0_v-0.0001) + m.integral(x1_v) + m.integral(x1_h) + m.integral(x0_h))

# m.options.OTOL = 1e-6
# m.options.RTOL = 1e-6

# m.options.IMODE = 6 # optimal control mode

m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True) # solve


virus_output = np.array(x0_v.value)
# host_output = np.array(x0_h.value)

# print(np.max(virus_output) / np.max(host_output))

plt.figure(1) # plot results
# plt.plot(m.time,x1.value,'k-',label=r'$x_1$')
# plt.plot(m.time,x2.value,'b-',label=r'$x_2$')
plt.plot(m.time*T0,virus_output,'b',label=r'$x0_v$')
# plt.plot(m.time*t_0,host_output,'r',label=r'$x0_h$')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Value')
plt.figure(2) # plot results
plt.plot(m.time,u0.value,'g',label=r'$u$')


plt.show()
