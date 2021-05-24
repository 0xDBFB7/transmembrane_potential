from transmembrane_lib import *
from scipy.integrate import odeint


host_cell = default_host_cell(np.array([]))
virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, np.array([]))



print(alpha_h)
print(beta_h)
print(gamma_h)
print(phi_h)
print(xi_h)


from sympy.abc import t # x is the independent variable
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp



u = Function('u')
t, alpha_v, beta_v, gamma_v, phi_v, xi_v = symbols("t alpha_v beta_v gamma_v phi_v xi_v")
alpha_h, beta_h, gamma_h, phi_h, xi_h = symbols("alpha_h beta_h gamma_h phi_h xi_h")



diffeq_one = Eq(xi_v*(1/(1+exp(-t))), (alpha_v*Derivative(u(t),t,t) + beta_v*Derivative(u(t),t) + gamma_v*u(t) - phi_v*(exp(-t))/((1+exp(-t))**2) - ((2*exp(-2*t))/((1+exp(-t))**3) - (exp(-t))/((1+exp(-t))**2))))
# diffeq_two = Eq(u(t), 0)

# system = [diffeq_one, diffeq_two]
dsolve(diffeq_one, u(t))


#
