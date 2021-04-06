from sympy import *
init_printing()
l_i, l_o, l_m, e_i, e_o, e_m, d, R, t  = symbols('l_i l_o l_m e_i e_o e_m d R t')

f = Function("f")
A_i = l_i + e_i * Derivative(f(t), t)
A_m = l_m + e_m * Derivative(f(t), t)
A_o = l_o + e_o * Derivative(f(t), t)

upper = 3 * A_o * (3 * d * R**2 * A_i + (2 * d**2 * R - d**3) * (A_m - A_i))
lower = 2 * R**3 * (A_m + 2 * A_o) * (A_m + (1/2) * A_i) - (2 * (R-d)**3)*(A_o - A_m)*(A_i-A_m)
#f = Function('f')
#pprint(dsolve(Derivative(f(t), t, t) + 9*f(t), f(t)))
#dsolve(Derivative(f(x), x, x) + 9*f(x), f(x))
#eq = (Eq(Derivative(x(t),t), 12*t*x(t) + 8*y(t)), Eq(Derivative(y(t),t), 21*x(t) + 7*t*y(t)))
#dsolve(eq)

dsolve(upper/lower,f(t))

#how about if we integrate in complex-freq space, then 
