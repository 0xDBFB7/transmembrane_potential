import sys
sys.path.append('../')
from transmembrane_lib import *
from scipy.integrate import odeint
from sympy.plotting import plot
from sympy import init_printing
import sympy
from sympy.abc import t
from sympy import Array, Sum, Indexed, IndexedBase, Idx

#export SYMPY_DEBUG=True

init_printing()
from sympy.abc import t # x is the independent variable
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp, pi

D = Derivative

u = Function('u')
P = Function('P')
lamda = Function('lamda')
t, alpha_v, beta_v, gamma_v, phi_v, xi_v,C1,C2,M = symbols("t alpha_v beta_v gamma_v phi_v xi_v C1 C2  M")
p_0, p_1, p_2, p_3, p_4, p_5 = symbols("p_0 p_1 p_2 p_3 p_4 p_5")
a_1, b_1 = symbols("a_1, b_1")
t_f = symbols("t_f")

j,m = symbols("j m", cls=Idx)
# alpha_h, beta_h, gamma_h, phi_h, xi_h = symbols("alpha_h beta_h gamma_h phi_h xi_h")


# p = Array([p_0, p_1, p_2, p_3, p_4, p_5])
# p_j = Function('p')(j)
# p = IndexedBase('p')
# P = Sum(p[j]*t**j,(j,0,5))#sympy chokes on all of these
P = p_0*t**0 + p_1*t**1 + p_2*t**2 + p_3*t**3 + p_4*t**4 + p_5*t**5
sympy.pprint(P)

# a_m = Function('a')(m)
# b_m = Function('b')(m)

a = IndexedBase('a')
b = IndexedBase('b')

# lamda = a_1 * cos(2 * 1 * pi * t / t_f) + b_1 * sin(2 * 1 * pi * t / t_f)
lamda = Sum(a[m] * cos(2 * 1 * pi * t / t_f), (m, 1, M)) + Sum(b[m] * sin(2 * 1 * pi * t / t_f), (m, 1, M))

# rhs = (P + lamda) #+ Derivative(P + lamda, t) + Derivative(P + lamda, t, t)
rhs = P + lamda

lhs = (alpha_v*Derivative(u(t),t,t) + beta_v*Derivative(u(t),t) + gamma_v*u(t))

solution = dsolve(Eq(rhs, lhs), u(t), ics={u(0): 0})

#u(t).diff(t).subs(t, 1e-6): 1
sympy.pprint(solution)

#conclusion: it is definitely analytic to go from one to the other.
