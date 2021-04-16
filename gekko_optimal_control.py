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

originally (mistakenly) used A6c (1998). Switched to eq 8, (multiplied by R, eq 10).
also missed switching the top and bottom: the top goes with the input terms (right) and the bottom with the output.

H(s)= (R*X(s)) / U(s) = (R a1 s^2 + R a2 s + R a3) / (b1 s^2 + b2 s + b3)

"Solution: Separate the equation so that the output terms, X(s), are on the
left and the input terms, Fa(s), are on the right.  Make sure there are only positive powers of s."

(b1 s^2 + b2 s + b3) X(s) = (R a1 s^2 + R a2 s + R a3) U(s)

"Now take the inverse Laplace Transform (so multiplications by "s" in the Laplace domain are replaced by derivatives in time)."

b1 x'' + b2 x' + b3 x = R a1 u'' + R a2 u' + R a3 u

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

    b1 x'' + b2 x' + b3 x = R a1 u'' + R a2 u' + R a3 u
    x''  = (R a1 u'' + R a2 u' + R a3 u - b2 x' - b3 x)/b1
    x2  = (R*a1*u2 + R*a2*u1 + R*a3*u0 - b2*x1 - b3*x0)/b1

docs: "In all simulation modes (IMODE=1,4,7), the number of equations must equal the number of variables."







"""


m = GEKKO() # initialize gekko
m.options.MAX_ITER = 100000
nt = 401
end = 1e-5
m.time = np.linspace(0,end,nt)
# Variables
x0_v = m.Var(value=0)
x1_v = m.Var()
x2_v = m.Var()

x0_h = m.Var(value=0)
x1_h = m.Var()
x2_h = m.Var()

t = m.Param(value=m.time)



u0 = m.Var()

# m.Equation(u0 == m.sin(t)) # for simulation
# m.Equation(u0 == m.exp(-(((t-(end/2.0))**2.0)/(2.0*((0.1e-4)**2.0))))) # for simulation


u1 = m.Var()
m.Equation(u1==u0.dt())
u2 = m.Var()
m.Equation(u2==u1.dt())


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

# m.Equation(x0_v == (R_v*a1_v*u2 + R_v*a2_v*u1 + R_v*a3_v*u0 - b2_v*x1_v  - x2_v*b1_v)/b3_v)
m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == (R_v*1e6*a1_v*u2 + R_v*1e6*a2_v*u1 + R_v*1e6*a3_v*u0 - b2_v*x1_v - b3_v*x0_v)/b1_v)
# #

m.Equation(x1_h==x0_h.dt())
m.Equation(x2_h==x1_h.dt())
m.Equation(x2_h == (R_h*1e6*a1_h*u2 + R_h*1e6*a2_h*u1 + R_h*1e6*a3_h*u0 - b2_h*x1_h - b3_h*x0_h)/b1_h)

#
# int_h = m.Var()
# m.Equation(int_h==m.integral(x0_h**2.0))
# m.Equation(t==m.vsum(m.abs2(u0)))

# m.Equation(m.integral(m.abs2(x0_h))==end)
# m.Equation(m.vsum(m.abs2(u0))==0.1)
# m.Equation(m.integral(m.abs2(u0))==0.1)
# integral()
# abs2()

m.Obj(-m.integral(m.abs2(x0_v))) # Objective function
m.options.IMODE = 6 # optimal control mode

# m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True) # solve


virus_output = np.array(x0_v.value)
host_output = np.array(x0_h.value)

print(virus_output)

print(np.max(virus_output) / np.max(host_output))


plt.figure(1) # plot results
# plt.plot(m.time,x1.value,'k-',label=r'$x_1$')
# plt.plot(m.time,x2.value,'b-',label=r'$x_2$')
plt.plot(m.time,virus_output,'b',label=r'$x0_v$')
plt.plot(m.time,host_output,'r',label=r'$x0_h$')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Value')
plt.figure(2) # plot results
plt.plot(m.time,u0.value,'g',label=r'$u$')



plt.show()
