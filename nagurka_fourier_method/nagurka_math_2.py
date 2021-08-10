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


t, 𝛂, 𝛃, 𝛄, 𝛙, 𝛏,C1,C2,M, J = symbols("t alpha beta gamma phi xi C1 C2  M J")
p_0, p_1, p_2, p_3, p_4, p_5 = symbols("p_0 p_1 p_2 p_3 p_4 p_5")
a_1, b_1 = symbols("a_1, b_1")
t_f = symbols("t_f")
j,m,k = symbols("j m k")

P = Function('P')
p = IndexedBase('p')
P = Sum(p[j]*t**j,(j,0,5))#sympy chokes on all of these


λ = Function('lamda')
a = IndexedBase('a')
b = IndexedBase('b')

λ = Sum(a[m] * cos(2 * m * pi * t / t_f), (m, 1, M)) + Sum(b[m] * sin(2 * m * pi * t / t_f), (m, 1, M))

sympy.pprint(λ)
sympy.pprint(diff(λ, t))
sympy.pprint(diff(λ, t, t))

X = P + λ

sympy.pprint(X )
sympy.pprint(diff(X , t))
sympy.pprint(diff(X , t, t))



'''


The fourier series and the polynomial define one time-course of one of the x outputs.

We want to obtain the u control input (and the other x output, of course), and then
the cost function J.
Unfortunately, both sympy and maxima choke on solving the DE for U after the fourier series,
even though this DE really should be analytic.

One could easily numerically integrate the problem, as discussed in one of Nagurka and Yen's examples
- but that's boring and kinda defeats the purpose of using this high-accuracy analytic method.



x2_v = 𝛂ᵥ*u2 + 𝛃ᵥ*u1 + 𝛄ᵥ*u0 - 𝛙ᵥ*x1_v - 𝛏ᵥ*x0_v




d²xᵥ/dt² = 𝛂*d²u/dt² + 𝛃*du/dt + 𝛄*u - 𝛙*dxᵥ/dt - 𝛏*xᵥ
𝛏*xᵥ = 𝛂*d²u/dt² + 𝛃*du/dt + 𝛄*u - 𝛙*dxᵥ/dt - d²xᵥ/dt²


𝛏xᵥ = 𝛂ü + 𝛃u̇ + 𝛄u - 𝛙ẋᵥ - ẍᵥ


We approximate xᵥ by Nagurka, Eq. 13:

xᵥ(t) = P(t) + λ(t)

    5
P = ∑ p_j
   j=0

    M
λ = ∑ a cos + b sin
   m=1

there is a subtlty here in that the derivative of a fourier series may not equal the term-by-term derivative.
https://math.stackexchange.com/questions/1754033/integration-and-differentiation-of-fourier-series
this is explicitly addressed in nagurka. this might be why the CASes choked - the terms that
ensure differentiability weren't known to it


𝛂ü + 𝛃u̇ + 𝛄u = 𝛏(P + λ)  + 𝛙ẋᵥ + ẍᵥ


now that's nothing more than the harmonic oscillator.

http://web.uvic.ca/~monahana/ode_1_flowchart.pdf
http://web.uvic.ca/~monahana/ode_2_flowchart.pdf
a brilliant ode flowchart!

this is a nonhomogenous differential equation.
Is the ODE linear? Yes.

𝛂ü + 𝛃u̇ + 𝛄u = 𝛏(P + λ) + 𝛙d(P + λ)/dt + d²(P + λ)/dt²

let us put it in standard form:

ü + 𝛃/𝛂u̇ + 𝛄/𝛂u = (𝛏/𝛂)(P + λ) + (𝛙/𝛂)d(P + λ)/dt + (1/𝛂)d²(P + λ)/dt²


The ODE is non-homogeneous.

"The solution is of form y(x) = yp(x) + yc(x) where yp(x) is a particular solution and yc(x) is
the solution of the associated homogeneous ODE."

"Find yc(x) by going back to Step 1." Homogeneous ODE:

ü + 𝛃/𝛂u̇ + 𝛄/𝛂u = 0

define two new constants,

A = 𝛃/𝛂
B = 𝛄/𝛂

typical values for
𝛂=1.3e-08
𝛃=5.7
𝛄=4e7

A² - 4B = 1e17
therefore, case A:

yc(t) = c₁e^(m₁t) + c₂e^(m₂t)


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
λU = Function('lamda')

λU = Sum(Aun[m] * a[m] * cos(2 * m * pi * t / t_f), (m, 1, M)) + Sum(Bun[m] * b[m] * sin(2 * m * pi * t / t_f), (m, 1, M))

U = PU + λU

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
#     RHS = sympy.simplify(diff(U,t,t) + (𝛃/𝛂) * diff(U,t) + (𝛄/𝛂)*U)
#     LHS = sympy.simplify((𝛏/𝛂)*X + (𝛙/𝛂)*diff(X,t) + (1/𝛂)*diff(X,t,t))
#
#     dill.dump_session(filename)


# Example 4,
# https://tutorial.math.lamar.edu/classes/de/undeterminedcoefficients.aspx
# "Now, as we’ve done in the previous examples we will need the coefficients of
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



So both sides are exactly the same barring the (𝛏/𝛂) and the undetermined coefficients.
maybe this is simpler than I thought (specifically because we're lucky enough to have the exact same form
on the input and output, just with different coefficients)

diff(U,t,t) + (𝛃/𝛂) * diff(U,t) + (𝛄/𝛂)*U = (𝛏/𝛂)*X + (𝛙/𝛂)*diff(X,t) + (1/𝛂)*diff(X,t,t)

do the undetermined coefficients have to be inside the sum (i.e. with an index?)
https://tutorial.math.lamar.edu/classes/calci/summationnotation.aspx
I think they can probably be outside.

(re-ordered RHS for convenience)
Q*diff(X,t,t) + (𝛃/𝛂)*R*diff(X,t) + (𝛄/𝛂)*S*X = (1/𝛂)*diff(X,t,t) + (𝛙/𝛂)*diff(X,t) + (𝛏/𝛂)*X
(...is that right?)

Q = (1/𝛂)

(𝛃/𝛂)*R = (𝛙/𝛂)
R = (𝛙/𝛂)/(𝛃/𝛂)

S = (𝛏/𝛂)/(𝛄/𝛂)


wait, but what is X?
no, this is all wrong.
'''






'''
integration:
https://math.stackexchange.com/questions/188567/when-can-a-series-be-integrated-term-by-term

'''



'''




₀₁₂₃₄₅

²³
cool unicode combining derivative,
combining third derivative
a bit small though

https://github.com/kragen/xcompose - maths-base
'''
