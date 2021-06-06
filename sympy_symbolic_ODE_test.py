from transmembrane_lib import *
from scipy.integrate import odeint


host_cell = default_host_cell(np.array([]))
virus = default_virus(np.array([]))

from sympy.abc import t
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp



u = Function('u')
alpha_v, beta_v, gamma_v, phi_v, xi_v = symbols("t alpha_v beta_v gamma_v phi_v xi_v")
x0_v = Function("x0_v")

alpha_h, beta_h, gamma_h, phi_h, xi_h = symbols("alpha_h beta_h gamma_h phi_h xi_h")
x0_h = Function("x0_h")



diffeq_one = Eq(Derivative(x0_v,t,t),
        alpha_v*Derivative(u(t),t,t) + beta_v*Derivative(u(t),t) + gamma_v*u(t) - phi_v*Derivative(x0_v,t) - xi_v*x0_v(t))
# diffeq_two = Eq(xi_v*(1/(1+exp(-t))), (*Derivative(u(t),t,t) + beta_v*Derivative(u(t),t) + gamma_v*u(t) - phi_v*(exp(-t))/((1+exp(-t))**2) - ((2*exp(-2*t))/((1+exp(-t))**3) - (exp(-t))/((1+exp(-t))**2))))

# diffeq_two = Eq(u(t), 0)

# system = [diffeq_one, diffeq_two]
dsolve(diffeq_one, u(t), ics={u(0): 0, x0_h(0): 0, x0_h(t).diff(t).subs(t, 0): 0, x0_v(0): 0, x0_v(t).diff(t).subs(t, 0): 0})


#
