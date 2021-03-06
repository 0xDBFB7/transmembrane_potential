'''
N.B.

Contains:

Derivatives of P and L; attempts at solving the whole diffeq using method of undetermined coefficients

'''
import dill
from pytexit import py2tex
import pickle

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
from sympy import Function, dsolve, Eq, Derivative, sin, cos, symbols, exp, pi, diff, Poly


t, 饾泜, 饾泝, 饾泟, 饾洐, 饾洀,C1,C2,M, J = symbols("t alpha beta gamma phi xi C1 C2  M J")
p_0, p_1, p_2, p_3, p_4, p_5 = symbols("p_0 p_1 p_2 p_3 p_4 p_5")
a_1, b_1 = symbols("a_1, b_1")
t_f = symbols("t_f")
j,m,k = symbols("j m k")

P = Function('P')
p = IndexedBase('p')
P = Sum(p[j]*t**j,(j,0,5))#sympy chokes on all of these


位 = Function('lamda')
a = IndexedBase('a')
b = IndexedBase('b')

位 = Sum(a[m] * cos(2 * m * pi * t / t_f), (m, 1, M)) + Sum(b[m] * sin(2 * m * pi * t / t_f), (m, 1, M))

sympy.pprint(位)
sympy.pprint(diff(位, t))
sympy.pprint(diff(位, t, t))

X = P + 位

sympy.pprint(X )
sympy.pprint(diff(X , t))
sympy.pprint(diff(X , t, t))



'''


In the inverse-dynamic use by Yen&Nagurka,
the fourier series and the polynomial define one time-course of one of the x outputs.

We want to obtain the u control input (and the other x output, of course), and then
the cost function J.
Unfortunately, both sympy and maxima choke on solving the DE for U after the fourier series,
even though this DE really should be analytic.

One could easily numerically integrate the problem, as discussed in one of Nagurka and Yen's examples
- but that's boring and kinda defeats the purpose of using this high-accuracy analytic method.



x2_v = 饾泜岬?*u2 + 饾泝岬?*u1 + 饾泟岬?*u0 - 饾洐岬?*x1_v - 饾洀岬?*x0_v




d虏x岬?/dt虏 = 饾泜*d虏u/dt虏 + 饾泝*du/dt + 饾泟*u - 饾洐*dx岬?/dt - 饾洀*x岬?
饾洀*x岬? = 饾泜*d虏u/dt虏 + 饾泝*du/dt + 饾泟*u - 饾洐*dx岬?/dt - d虏x岬?/dt虏


饾洀x岬? = 饾泜u虉 + 饾泝u虈 + 饾泟u - 饾洐x虈岬? - x虉岬?


We approximate x岬? by Nagurka, Eq. 13:

x岬?(t) = P(t) + 位(t)

    5
P = 鈭? p_j
   j=0

    M
位 = 鈭? a cos + b sin
   m=1

there is a subtlty here in that the derivative of a fourier series may not equal the term-by-term derivative.
https://math.stackexchange.com/questions/1754033/integration-and-differentiation-of-fourier-series
this is explicitly addressed in nagurka. this might be why the CASes choked - the terms that
ensure differentiability weren't known to it


饾泜u虉 + 饾泝u虈 + 饾泟u = 饾洀(P + 位)  + 饾洐x虈岬? + x虉岬?


now that's nothing more than the harmonic oscillator.

http://web.uvic.ca/~monahana/ode_1_flowchart.pdf
http://web.uvic.ca/~monahana/ode_2_flowchart.pdf
a brilliant ode flowchart!

this is a nonhomogenous differential equation.
Is the ODE linear? Yes.

饾泜u虉 + 饾泝u虈 + 饾泟u = 饾洀(P + 位) + 饾洐d(P + 位)/dt + d虏(P + 位)/dt虏

let us put it in standard form:

u虉 + 饾泝/饾泜u虈 + 饾泟/饾泜u = (饾洀/饾泜)(P + 位) + (饾洐/饾泜)d(P + 位)/dt + (1/饾泜)d虏(P + 位)/dt虏


The ODE is non-homogeneous.

"The solution is of form y(x) = yp(x) + yc(x) where yp(x) is a particular solution and yc(x) is
the solution of the associated homogeneous ODE."

"Find yc(x) by going back to Step 1." Homogeneous ODE:

u虉 + 饾泝/饾泜u虈 + 饾泟/饾泜u = 0

define two new constants,

A = 饾泝/饾泜
B = 饾泟/饾泜

typical values for
饾泜=1.3e-08
饾泝=5.7
饾泟=4e7

A虏 - 4B = 1e17
therefore, case A:

yc(t) = c鈧乪^(m鈧乼) + c鈧俥^(m鈧倀)


(I knew I would mess up beta and B, man, lack of CAS is hurting).

(spoiler alert: We actually didn't need yc at all, it seems)

"Is g(x) a sum/product of constants, polynomials, exponentials, and sine/cosine
functions?" Yes.
"Can find yp(x) using method of undetermined coefficients. Let the functions hi(x) be made
up of the functions in g(x) and their derivatives.
Do any of the hi(x) belong to the fundamental set y1(x), y2(x) of the homogeneous ODE?"

https://tutorial.math.lamar.edu/classes/de/undeterminedcoefficients.aspx


'''


# our guess at a solution with undetermined coefficients
Aun = IndexedBase('Aun') #undetermined
Bun = IndexedBase('Bun')
Cun = IndexedBase('Cun')

PU = Function('PU')
PU = Sum(Cun[j]*p[j]*t**j,(j,0,5))
位U = Function('lamda')

位U = Sum(Aun[m] * a[m] * cos(2 * m * pi * t / t_f), (m, 1, M)) + Sum(Bun[m] * b[m] * sin(2 * m * pi * t / t_f), (m, 1, M))

U = PU + 位U

# sympy.pprint(U )
# sympy.pprint(diff(U , t))
sympy.pprint(diff(X , t, t))
# sympy.pprint(diff(U , t, t))

print(sympy.latex(X))
print(sympy.latex(diff(X, (t))))
print(sympy.latex(diff(X , t, t)))

z = symbols("z",  integer = True)
print()
print(sympy.latex(diff(P, (t, z)).doit()))
# doesn't work - sympy doesn't know how to take this arbitrary deterivative?


print(sympy.latex(diff(P, t, t, t, t, t, t)))


# print(sympy.latex(diff(U , t, t)))
# print()
# print(X)
# print(diff(X , t))
# print(diff(X , t, t))

# sympy.pprint(sympy.integrate(X, t))

#
# filename = 'cache/cache.pkl'
# try:
#     dill.load_session(filename)
# except:
#     RHS = sympy.simplify(diff(U,t,t) + (饾泝/饾泜) * diff(U,t) + (饾泟/饾泜)*U)
#     LHS = sympy.simplify((饾洀/饾泜)*X + (饾洐/饾泜)*diff(X,t) + (1/饾泜)*diff(X,t,t))
#
#     dill.dump_session(filename)


# Example 4,
# https://tutorial.math.lamar.edu/classes/de/undeterminedcoefficients.aspx
# "Now, as we鈥檝e done in the previous examples we will need the coefficients of
# the terms on both sides of the equal sign to be the same so set coefficients equal and solve."

# sympy.pprint(RHS)
# sympy.pprint(LHS)

# print(sympy.latex(Eq(RHS,LHS), mode="equation"))
# sympy.pprint(RHS.args[0])
# sympy.pprint(LHS.args[0])

#https://stackoverflow.com/questions/26927481/finding-coefficient-of-a-specific-term-using-sympy

# print(RHS.coeff(t))
# print(LHS.coeff(t))

# print(Poly(RHS, t).all_coeffs())

'''



So both sides are exactly the same barring the (饾洀/饾泜) and the undetermined coefficients.
maybe this is simpler than I thought (specifically because we're lucky enough to have the exact same form
on the input and output, just with different coefficients)

diff(U,t,t) + (饾泝/饾泜) * diff(U,t) + (饾泟/饾泜)*U = (饾洀/饾泜)*X + (饾洐/饾泜)*diff(X,t) + (1/饾泜)*diff(X,t,t)

do the undetermined coefficients have to be inside the sum (i.e. with an index?)
https://tutorial.math.lamar.edu/classes/calci/summationnotation.aspx
I think they can probably be outside.


(re-ordered RHS for convenience)
Q*diff(X,t,t) + (饾泝/饾泜)*R*diff(X,t) + (饾泟/饾泜)*S*X = (1/饾泜)*diff(X,t,t) + (饾洐/饾泜)*diff(X,t) + (饾洀/饾泜)*X
(...is that right?)

Q = (1/饾泜)

(饾泝/饾泜)*R = (饾洐/饾泜)
R = (饾洐/饾泜)/(饾泝/饾泜)

S = (饾洀/饾泜)/(饾泟/饾泜)

wait, but what is X?
no, this is all wrong.
'''


"""
Much later:

Introducing two new variables R and W as placeholders for U and X, whichever order is chosen,


SymPy doesn't like floating-point exponents! doesn't work at all!
"""

a1,b1,a2,b2,p_0,p_1,p_2,p_3,p_4,p_5 = symbols("a1 b1 a2 b2 p_0 p_1 p_2 p_3 p_4 p_5", real=True)
R = cos(a1*t) #p_2*t**2 #cos(a1*t) #+ sin(b1*t) #cos(a2*t) + sin(b2*t) + p0 + p1*t
# U = sin(t) + cos(t)
W = Function('W') #previously U
A,B,C,D,E = symbols("A,B,C,D,E",real=True)
LHS = diff(W(t),t,t) + A * diff(W(t),t) + B*W(t)
RHS = C*diff(R,t,t) + D*diff(R,t) + E*R

sympy.pprint(Eq(LHS, RHS))

W_0, d_W_0, W_tf, d_W_tf = symbols("W_0 d_W_0 W_tf d_W_tf", real=True)

solution = dsolve(Eq(LHS, RHS), W(t), ics={W(0): 0, W(t).diff(t,1).subs(t,0): 0})

# adding W(t_f): W_tf, W(t).diff(t,1).subs(t,t_f): d_W_tf}
# raises ValueError: Couldn't solve for initial conditions
sympy.pprint(solution)
# sympy.pprint(diff(solution.args[1], t,t))
# sympy.pprint(sympy.simplify(solution))
print(len(solution.args[1].args)) #first arg is just u
# a1b1 single term, 18 + 2 for the IC particular
# a1b1+a2b2; 36 + 2

# sympy.pprint(sympy.integrate(solution, t))

# W = sympy.solve(solution, W(t))[0]



"""
So everything is nicely analytic, but the whole equation is algebraicly tricky.

takes maybe 20 minutes or something 
"""



"""
That produces a huge amount of algebra. Can the Laplace transform cut down on that?
http://josephcslater.github.io/solve-ode.html
"""

print("ilt----------------")
# s = sympy.symbols('s')
#
# LHS = W(s)*s**2 + A * W(s)*s + B*W(s)
# lt = Eq(LHS, sympy.laplace_transform(RHS, t, s, noconds = True))
# sympy.pprint(lt)
#
# Ult = sympy.solve(lt,W(s))
# sympy.pprint(Ult)
# sympy.pprint(sympy.inverse_laplace_transform(Ult[0],s,t))

"""
Nope! That's way less manageable. okay!
"""

"""
Integrating the pore N function is somewhat expensive.
can that be done analytically?
"""
# v_m,alpha,N0,q,v_ep, N_ic = symbols("v_m,alpha,N0,q,v_ep, N_ic", real=True)
# N = Function('N')
# # v_m = Function('v_m')

# v_m = sympy.simplify(solution.rhs)

# k = (v_m/v_ep)**2.0

# LHS = diff(N(t),t)
# RHS = alpha * exp(k) * (1.0 - (N(t)/N0)*exp(-q*k))

# sympy.pprint(Eq(LHS, RHS))

# solution = dsolve(Eq(LHS, RHS), N(t), ics={N(0): N_ic})

# # adding W(t_f): W_tf, W(t).diff(t,1).subs(t,t_f): d_W_tf}
# # raises ValueError: Couldn't solve for initial conditions
# sympy.pprint(solution)
# print(solution)
# # sympy.pprint(sympy.simplify(solution))


"""
Oh yeah, it can.


Hey, if everything's analytic and the cost function is just N, won't gradient methods work great?


"""

"""
we have the pore current equation (which doesn't depend on t explicitly).
dV/dt = i/c but we need the second derivative of V, so derive in terms of unspecified V(t)


"""

w0,F,R,T,r_m,diameter,h,sigma,n,N = symbols("w0,F,R,T,r_m,diameter,h,sigma,n,N")
v_m = Function('v_m')

# v_m = (V_m(t)) * (F/(R*T))
v_m = v_m(t)

i_ep_term_1 = (pi * (r_m**2) * sigma * v_m * R * T) / (F * h)

i_ep_term_2_divisor_1 = exp(v_m)*((w0*exp(w0-n*v_m) - (n*v_m)) / (w0 - (n*v_m)))
i_ep_term_2_divisor_2 = -((w0*exp(w0+n*v_m) + (n*v_m)) / (w0 + (n*v_m)))
i_ep_term_2 = exp(v_m - 1) / (i_ep_term_2_divisor_1 + i_ep_term_2_divisor_2)

i_ep = i_ep_term_1 * i_ep_term_2
# I = C_m dv/dt
# dv_dt = I/C_m
# permittivity already has the eps0 in it 
A = 4*pi*(diameter/2)**2
C_m = k * A / h
I_ep = -(i_ep * N / C_m)

sympy.pprint(I_ep)
sympy.pprint(I_ep.diff(t))
print(I_ep.diff(t))
# sympy.pprint(sympy.simplify(I_ep.diff(t)))
# this line takes too long. put it in SageMath, simplified with the maxima backend.
# that

# -N*r_m^2*sigma*((F*n*w0*e^(w0 - F*n*V_m(t)/(R*T))*diff(V_m(t), t)/(R*T) + F*n*
# diff(V_m(t), t)/(R*T))*e^(F*V_m(t)/(R*T))/(w0 - F*n*V_m(t)/(R*T)) - (w0*e^(w0 
# - F*n*V_m(t)/(R*T)) - F*n*V_m(t)/(R*T))*F*n*e^(F*V_m(t)/(R*T))*diff(V_m(t), t)
# /(R*T*(w0 - F*n*V_m(t)/(R*T))^2) - (w0*e^(w0 - F*n*V_m(t)/(R*T)) - F*n*V_m(t)/
# (R*T))*F*e^(F*V_m(t)/(R*T))*diff(V_m(t), t)/(R*T*(w0 - F*n*V_m(t)/(R*T))) + (F
# *n*w0*e^(w0 + F*n*V_m(t)/(R*T))*diff(V_m(t), t)/(R*T) + F*n*diff(V_m(t), t)/(R
# *T))/(w0 + F*n*V_m(t)/(R*T)) - (w0*e^(w0 + F*n*V_m(t)/(R*T)) + F*n*V_m(t)/(R*T
# ))*F*n*diff(V_m(t), t)/(R*T*(w0 + F*n*V_m(t)/(R*T))^2))*V_m(t)*e^(F*V_m(t)/(R*
# T) - 1)/(Cell_D^2*k*((w0*e^(w0 - F*n*V_m(t)/(R*T)) - F*n*V_m(t)/(R*T))*e^(F*V_
# m(t)/(R*T))/(w0 - F*n*V_m(t)/(R*T)) - (w0*e^(w0 + F*n*V_m(t)/(R*T)) + F*n*V_m(
# t)/(R*T))/(w0 + F*n*V_m(t)/(R*T)))^2) - N*r_m^2*sigma*e^(F*V_m(t)/(R*T) - 1)*d
# iff(V_m(t), t)/(Cell_D^2*k*((w0*e^(w0 - F*n*V_m(t)/(R*T)) - F*n*V_m(t)/(R*T))*
# e^(F*V_m(t)/(R*T))/(w0 - F*n*V_m(t)/(R*T)) - (w0*e^(w0 + F*n*V_m(t)/(R*T)) + F
# *n*V_m(t)/(R*T))/(w0 + F*n*V_m(t)/(R*T)))) - F*N*r_m^2*sigma*V_m(t)*e^(F*V_m(t
# )/(R*T) - 1)*diff(V_m(t), t)/(Cell_D^2*R*T*k*((w0*e^(w0 - F*n*V_m(t)/(R*T)) - 
# F*n*V_m(t)/(R*T))*e^(F*V_m(t)/(R*T))/(w0 - F*n*V_m(t)/(R*T)) - (w0*e^(w0 + F*n
# *V_m(t)/(R*T)) + F*n*V_m(t)/(R*T))/(w0 + F*n*V_m(t)/(R*T))))

# simplified by taking out the V_m -> v_m reduction:




'''
integration:
https://math.stackexchange.com/questions/188567/when-can-a-series-be-integrated-term-by-term

'''



'''




鈧?鈧佲倐鈧冣倓鈧?

虏鲁
cool unicode combining derivative,
combining third derivative
a bit small though

https://github.com/kragen/xcompose - maths-base
'''
