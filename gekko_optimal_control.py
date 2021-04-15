from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt

from transmembrane_lib import *

t = np.array([0])

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 50e-9, 5e-9, t)

virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, t)


"""

PSOPT_Manual_R5.pdf

For convenience, let's convert the Kotnik equation to standard DE form.

per https://lpsa.swarthmore.edu/Representations/SysRepTransformations/TF2SDE.html

originally (mistakenly) used A6c (1998):

H(s)/R= X(s) / F(s) = (R*a1 + R*a2 s + R*a3 s^2) / (b1 s + b2 s^2 + b3 s^3)

"Solution: Separate the equation so that the output terms, X(s), are on the
left and the input terms, Fa(s), are on the right.  Make sure there are only positive powers of s."

(R*a1 + R*a2 s + R*a3 s^2) X(s) = (b1 s + b2 s^2 + b3 s^3) U(s)

"Now take the inverse Laplace Transform (so multiplications by "s" in the Laplace domain are replaced by derivatives in time)."

(R*a1 x + R*a2 x' + R*a3 x'') = (b1 u' + b2 u'' + b3 u''')

where u is the input function,

now, gekko wants it in first order, so we have to substitute (Convert to a System of DEs)
http://www.math.utah.edu/~gustafso/2250systems-de.pdf
http://www.sharetechnote.com/html/DE_HigherOrderDEtoFirstOrderDE.html
https://math.berkeley.edu/~zworski/128/psol12.pdf

https://math.stackexchange.com/questions/1120984/

(reworked because I found the notation confusing (u3 should be the third derivative!))

u0 = u
u1 = u'
u2 = u''
u3 = u'''

      x0 = x
x0' = x1 = x'
x1' = x2 = x'' =

    R*a1 x + R*a2 x' + R*a3 x'' = (b1 u' + b2 u'' + b3 u''')
    x'' = (b1 u' + b2 u'' + b3 u''' - R*a1 x - R*a2 x') / R*a3
    x2 = (b1 u1 + b2 u2 + b3 u3 - R*a1 x0 - R*a2 x1) / R*a3

    x  = ((b1 u' + b2 u'' + b3 u''') - R*a2 x' - R*a3 x'') / R*a1

docs: "In all simulation modes (IMODE=1,4,7), the number of equations must equal the number of variables."







"""


m = GEKKO() # initialize gekko
nt = 101
m.time = np.linspace(0,1e-6,nt)
# Variables
x0_v = m.Var()
x1_v = m.Var()
x2_v = m.Var()

x0_h = m.Var(value=1)
x1_h = m.Var()
x2_h = m.Var()

t = m.Param(value=m.time)

int_h = m.Var()


u0 = m.Var()
m.Equation(u0 == m.sin(t)) # for simulation
u1 = m.Var()
m.Equation(u1==u0.dt())
u2 = m.Var()
m.Equation(u2==u1.dt())
u3 = m.Var()
m.Equation(u3==u2.dt())

# p = np.zeros(nt) # mark final time point
# p[-1] = 1.0
#
a1_v = m.Const(virus.a_1)
a2_v = m.Const(virus.a_2)
a3_v = m.Const(virus.a_3)
b1_v = m.Const(virus.b_1)
b2_v = m.Const(virus.b_2)
b3_v = m.Const(virus.b_3)
R_v = m.Const(virus.R)


a1_h = m.Const(host_cell.a_1)
a2_h = m.Const(host_cell.a_2)
a3_h = m.Const(host_cell.a_3)
b1_h = m.Const(host_cell.b_1)
b2_h = m.Const(host_cell.b_2)
b3_h = m.Const(host_cell.b_3)
R_h = m.Const(host_cell.R)


# Equations

# x2_v==((b1_v*u1 + b2_v*u2 + b3_v*u3 - R_v*a1_v*x0_v + R_v*a2_v*x1_v) / R_v*a3_v)
m.Equation(x0_v==(((b1_v*u1 + b2_v*u2 + b3_v*u3) - R_v*a2_v*x1_v - R_v*a3_v*x2_v) / R_v*a1_v))
m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==((b1_v*u1 + b2_v*u2 + b3_v*u3 - R_v*a1_v*x0_v - R_v*a2_v*x1_v) / R_v*a3_v))
# #

# m.Equation(x0_h==(((b1_h*u1 + b2_h*u2 + b3_h*u3) - R_h*a2_h*x1_h - R_h*a3_h*x2_h) / R_h*a1_h))
# m.Equation(x1_h==x0_h.dt())
# m.Equation(x2_h==((b1_h*u1 + b2_h*u2 + b3_h*u3 - R_h*a1_h*x0_h - R_h*a2_h*x1_h) / R_h*a3_h))

# m.Equation(int_h==m.integral(x0_h))


# integral()
# abs2()

m.Obj(-x0_v + int_h) # Objective function
m.options.IMODE = 6 # optimal control mode

# m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True) # solve
plt.figure(1) # plot results
# plt.plot(m.time,x1.value,'k-',label=r'$x_1$')
# plt.plot(m.time,x2.value,'b-',label=r'$x_2$')
plt.plot(m.time,u0.value,'r--',label=r'$u$')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()
