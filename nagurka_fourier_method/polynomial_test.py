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

# To obtain an expression for the boundary conditions of a 5-term polynomial

t, t_f = symbols("t t_f")
P, p_0, p_1, p_2, p_3, p_4, p_5 = symbols("P p_0 p_1 p_2 p_3 p_4 p_5")
# P = Function('P')(t)
P_0, P_tf, d_P_0, d_P_tf, d_d_P_0, d_d_P_tf = symbols("P_0, P_tf, d_P_0, d_P_tf, d_d_P_0, d_d_P_tf ")

P = p_0 + p_1*t + p_2*t**2 + p_3*t**3 + p_4*t**4 + p_5*t**5

Exp = [ (P - P_0).subs(t, 0),
        (P - P_tf).subs(t, t_f),
        (diff(P, t) - d_P_0).subs(t, 0),
        (diff(P, t) - d_P_tf).subs(t, t_f),
        (diff(P, t, t) - d_d_P_0).subs(t, 0),
        (diff(P, t, t) - d_d_P_tf).subs(t, t_f)]

sympy.pprint(Exp)
sympy.pprint(sympy.solve(Exp, p_0, p_1, p_2, p_3, p_4, p_5))
                                # Eq(P(0), P0), Eq(P(t_f), P0)]))
