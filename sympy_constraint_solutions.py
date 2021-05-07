from transmembrane_lib import *
from scipy.integrate import odeint
from sympy import init_printing
import sympy

init_printing()

'''
What is the family of solutions for which x_h is zero?
'''

#
# print(alpha_h)
# print(beta_h)
# print(gamma_h)
# print(phi_h)
# print(xi_h)


from sympy.abc import t # x is the independent variable
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp



u = Function('u')
t, alpha_v, beta_v, gamma_v, phi_v, xi_v = symbols("t alpha_v beta_v gamma_v phi_v xi_v")
alpha_h, beta_h, gamma_h, phi_h, xi_h = symbols("alpha_h beta_h gamma_h phi_h xi_h")

solution = dsolve(Eq(0, (alpha_h*Derivative(u(t),t,t) + beta_h*Derivative(u(t),t) + gamma_h*u(t))), u(t),ics={u(0): 0})

sympy.pprint(solution)
print(str(solution))

#
