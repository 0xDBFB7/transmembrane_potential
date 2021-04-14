from gekko import GEKKO
import numpy as np
import matplotlib.pyplot as plt




"""

PSOPT_Manual_R5.pdf

For convenience, let's convert the Kotnik equation to standard DE form.

per https://lpsa.swarthmore.edu/Representations/SysRepTransformations/TF2SDE.html

A6c (1998):

H(s)/R= X(s) / F(s) = (R*a1 + R*a2 s + R*a3 s^2) / (b1 s + b2 s^2 + b3 s^3)

"Solution: Separate the equation so that the output terms, X(s), are on the
left and the input terms, Fa(s), are on the right.  Make sure there are only positive powers of s."

(R*a1 + R*a2 s + R*a3 s^2) X(s) = (b1 s + b2 s^2 + b3 s^3) U(s)

"Now take the inverse Laplace Transform (so multiplications by "s" in the Laplace domain are replaced by derivatives in time)."

(R*a1 x + R*a2 x' + R*a3 x'') = (b1 u' + b2 u'' + b3 u''')

where u is the input function,

now, gekko wants it in first order, so we have to substitute:
http://www.sharetechnote.com/html/DE_HigherOrderDEtoFirstOrderDE.html
https://math.stackexchange.com/questions/1120984/

https://math.berkeley.edu/~zworski/128/psol12.pdf



"""



m = GEKKO() # initialize gekko
nt = 101
m.time = np.linspace(0,2,nt)
# Variables
x1 = m.Var(value=1)
x2 = m.Var(value=0)
u = m.Var(value=0,lb=-1,ub=1)
p = np.zeros(nt) # mark final time point
p[-1] = 1.0
final = m.Param(value=p)
# Equations
m.Equation(x1.dt()==u)
m.Equation(x2.dt()==0.5*x1**2)

integral()
abs2()


m.Obj(x2*final) # Objective function
m.options.IMODE = 6 # optimal control mode
m.solve(disp=False) # solve
plt.figure(1) # plot results
plt.plot(m.time,x1.value,'k-',label=r'$x_1$')
plt.plot(m.time,x2.value,'b-',label=r'$x_2$')
plt.plot(m.time,u.value,'r--',label=r'$u$')
plt.legend(loc='best')
plt.xlabel('Time')
plt.ylabel('Value')
plt.show()
