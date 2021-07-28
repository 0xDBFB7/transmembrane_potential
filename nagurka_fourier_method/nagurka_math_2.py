import sys
sys.path.append('../')
from transmembrane_lib import *
from scipy.integrate import odeint
from sympy.plotting import plot
from sympy import init_printing
import sympy
from sympy.abc import t
from sympy import Array, Sum, Indexed, IndexedBase, Idx
init_printing()
from sympy.abc import t # x is the independent variable
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp, pi, diff



t, alpha_v, beta_v, gamma_v, phi_v, xi_v,C1,C2,M = symbols("t alpha_v beta_v gamma_v phi_v xi_v C1 C2  M")
p_0, p_1, p_2, p_3, p_4, p_5 = symbols("p_0 p_1 p_2 p_3 p_4 p_5")
a_1, b_1 = symbols("a_1, b_1")
t_f = symbols("t_f")
j,m,k = symbols("j m k", cls=Idx)

lamda = Function('lamda')
a = IndexedBase('a')
b = IndexedBase('b')

lamda = Sum(a[m] * cos(2 * m * pi * t / t_f), (m, 1, M)) + Sum(b[m] * sin(2 * m * pi * t / t_f), (m, 1, M))

sympy.pprint(lamda)
sympy.pprint(diff(lamda, t))
sympy.pprint(diff(lamda, t, t))
print(lamda)
print(diff(lamda, t))
print(diff(lamda, t, t))
