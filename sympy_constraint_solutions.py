from transmembrane_lib import *
from scipy.integrate import odeint
from sympy.plotting import plot
from sympy import init_printing
import sympy

init_printing()

'''
What is the family of solutions for which x_h is zero?
'''

host_cell = default_host_cell(np.zeros(0))

from sympy.abc import t # x is the independent variable
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp



u = Function('u')
t, alpha_v, beta_v, gamma_v, phi_v, xi_v,C1,C2 = symbols("t alpha_v beta_v gamma_v phi_v xi_v, C1, C2")
# alpha_h, beta_h, gamma_h, phi_h, xi_h = symbols("alpha_h beta_h gamma_h phi_h xi_h")


solution = dsolve(Eq(0, (host_cell.alpha*Derivative(u(t),t,t) + host_cell.beta*Derivative(u(t),t) + host_cell.gamma*u(t))), u(t),
                ics={u(0): 0})

#u(t).diff(t).subs(t, 1e-6): 1
sympy.pprint(solution)

solution = solution.rhs
solution_first_derivative = sympy.diff(solution, t)
solution_second_derivative = sympy.diff(solution_first_derivative,t)

sum_unsub = solution + solution_first_derivative + solution_second_derivative

sum = solution.subs(t, 1e-9) + solution_first_derivative.subs(t, 1e-9) + solution_second_derivative.subs(t, 1e-9)

sympy.pprint(sum.subs(C2,0))
plot(sum.subs(C2,0),(t,0,1e-6))

'''
Conclusion:

the solution produced only equals 0 when C1/C2=0. (as confirmed by _analytic_constraint),
which is useless.

'''
