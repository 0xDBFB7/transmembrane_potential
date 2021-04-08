from dataclasses import dataclass
from math import pi, sqrt, e, log, isclose, exp
# from scipy.optimize import curve_fit
import numpy as np
from numpy import heaviside
# import matplotlib.pyplot as plt
# import h5py
from scipy.constants import epsilon_0, mu_0
from pytexit import py2tex

def ustep(v):
    # v[v>0] = v[v>0]
    # v[v<0] = 0
    # return v
    # if(v > 0.0):
    #     return v
    # else:
    #     return 0.0
    return heaviside(v, 0.5)

@dataclass
class Cell:
    extracellular_conductivity: float # S/m
    extracellular_permittivity: float # relative
    intracellular_conductivity: float # S/m
    intracellular_permittivity: float # relative
    membrane_conductivity: float # S/m
    membrane_permittivity: float # relative
    cell_diameter: float # meters
    membrane_thickness: float

'''
This is a verbatim implementation of
[1]:
Kotnik T, Miklavčič D, Slivnik T.
Time course of transmembrane voltage induced by time-varying electric fields—a method for theoretical analysis and its application.
Bioelectrochemistry and Bioenergetics 1998;45:3–16.
https://doi.org/10.1016/S0302-4598(97)00093-7.

Kotnik et al also have number of other papers with minor variations that may be of some use

If a truly arbitrary
the equivalent of the discrete fourier transform for the Laplace transform appears to be the Z-transform.

'''

def tau_1_f(b_1, b_2, b_3):
    #equation A6e in Kotnik
    return (2.0 * b_3) / (b_2 - sqrt(((b_2)**2.0) - ((4.0*b_1) * b_3)))

def tau_2_f(b_1, b_2, b_3):
    #equation A6f in Kotnik
    return (2.0 * b_3) / (b_2 + sqrt(((b_2)**2.0) - ((4.0*b_1) * b_3)))



def delta_transmembrane_unit_ramp(t, cell):

    e_o = cell.extracellular_permittivity * epsilon_0 # S/m
    e_i = cell.intracellular_permittivity * epsilon_0 #S/m
    e_m = cell.membrane_permittivity * epsilon_0 #S/m
    R = cell.cell_diameter / 2.0

    l_o = cell.extracellular_conductivity # S/m
    l_i = cell.intracellular_conductivity #S/m
    l_m = cell.membrane_conductivity #S/m

    d = cell.membrane_thickness
    # epsilon_0

    sub1 = (3.0 * (R**2.0) - 3.0 * d * R + d**2.0)
    sub2 = (3.0 * d * R - d**2.0)

    a_1 = 3.0 * d * l_o * ((l_i * (sub1)) + l_m*(sub2)) #eq.9a
    a_2 = 3.0 * d * ((l_i * e_o + l_o * e_i) * sub1 + (l_m * e_o + l_o * e_m) * sub2)
    a_3 = 3.0 * d * e_o * (e_i * (sub1) + e_m * sub2)

    b_1 = 2.0 * R**3.0 * (l_m +     2.0*l_o) * (l_m + 0.5 * l_i) + 2.0 * (R-d)**3.0 * (l_m - l_o) * (l_i - l_m)

    b_2 = 2.0 * R**3.0 * (l_i * (0.5 * e_m + e_o) + l_m * (0.5*e_i + 2.0*e_m + 2*e_o) + l_o * (e_i + 2.0 * e_m)) + (2.0 * (R - d)**3.0\
    * (l_i * (e_m - e_o) + l_m * (e_i - 2.0*e_m + e_o) - l_o * (e_i - e_m))) # is this truly a multiply, or a cross?


    b_3 = 2.0 * R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)

    #py2tex("b_3 = 2.0 - R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)")
    #https://quicklatex.com/

    # Kotnik variously use "step function" or the "unit step".

    tau_1 = tau_1_f(b_1, b_2, b_3)
    tau_2 = tau_2_f(b_1, b_2, b_3)

    # a9d, Kotnik 1998, unit ramp function response
    delta_phi_m_t = (a_1 / b_1) * t * ustep(t)
    phisub1 = (a_2 / (2 * b_1)) - ((a_1*b_2) / (2 * (b_1 **2.0)))
    phisub2 = (((a_1 * b_3)/b_1) + ((a_2 * b_2) / (2.0 * b_1)) - ((a_1 * (b_2**2.0)) / (2.0 * (b_1**2.0))) - a_3)
    phisub3 = np.sqrt(b_2**2.0 - 4.0 * b_1 * b_3)

    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/tau_1)) * ustep(t)
    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/tau_2)) * ustep(t) # glitch?
    delta_phi_m_t *= R
    # # unlike in symbolic math, the unit step doesn't actually prevent errors from occuring in the np.exp step.
    delta_phi_m_t = np.nan_to_num(delta_phi_m_t)
    # if(np.isscalar(t)):
    #     if(t < 0):
    #         delta_phi_m_t == 0

    # if(np.isscalar(t)):


    return delta_phi_m_t

def delta_transmembrane_unit_step(t, cell):
    e_o = cell.extracellular_permittivity * epsilon_0 # S/m
    e_i = cell.intracellular_permittivity * epsilon_0 #S/m
    e_m = cell.membrane_permittivity * epsilon_0 #S/m
    R = cell.cell_diameter / 2.0

    l_o = cell.extracellular_conductivity # S/m
    l_i = cell.intracellular_conductivity #S/m
    l_m = cell.membrane_conductivity #S/m

    d = cell.membrane_thickness
    # epsilon_0

    sub1 = (3.0 * (R**2.0) - 3.0 * d * R + d**2.0)
    sub2 = (3.0 * d * R - d**2.0)

    a_1 = 3.0 * d * l_o * ((l_i * (sub1)) + l_m*(sub2)) #eq.9a
    a_2 = 3.0 * d * ((l_i * e_o + l_o * e_i) * sub1 + (l_m * e_o + l_o * e_m) * sub2)
    a_3 = 3.0 * d * e_o * (e_i * (sub1) + e_m * sub2)

    b_1 = 2.0 * R**3.0 * (l_m +     2.0*l_o) * (l_m + 0.5 * l_i) + 2.0 * (R-d)**3.0 * (l_m - l_o) * (l_i - l_m)

    b_2 = 2.0 * R**3.0 * (l_i * (0.5 * e_m + e_o) + l_m * (0.5*e_i + 2.0*e_m + 2*e_o) + l_o * (e_i + 2.0 * e_m)) + (2.0 * (R - d)**3.0\
    * (l_i * (e_m - e_o) + l_m * (e_i - 2.0*e_m + e_o) - l_o * (e_i - e_m))) # is this truly a multiply, or a cross?


    b_3 = 2.0 * R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)

    #py2tex("b_3 = 2.0 - R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)")
    #https://quicklatex.com/

    # Kotnik variously use "step function" or the "unit step".

    tau_1 = tau_1_f(b_1, b_2, b_3)
    tau_2 = tau_2_f(b_1, b_2, b_3)

    # a9d, Kotnik 1998, unit step function response
    delta_phi_m_t = (a_3 / b_3) * ustep(t)
    phisub1 = (a_1 / (2.0 * b_1)) - (a_3 / (2.0 * b_3))
    phisub2 = ((a_1 * b_2) /  (2.0 * b_1)) - a_2 + ((a_3 * b_2) / (2.0 * b_3))
    phisub3 = np.sqrt(b_2**2.0 - 4.0 * b_1 * b_3)

    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/tau_1)) * ustep(t)
    delta_phi_m_t += (phisub1 - (phisub2/phisub3)) * (1.0 - np.exp(-t/tau_2)) * ustep(t) # glitch?
    delta_phi_m_t *= R
    # # unlike in symbolic math, the unit step doesn't actually prevent errors from occuring in the np.exp step.
    delta_phi_m_t = np.nan_to_num(delta_phi_m_t)
    # if(np.isscalar(t)):
    #     if(t < 0):
    #         delta_phi_m_t == 0

    # if(np.isscalar(t)):


    return delta_phi_m_t


def delta_transmembrane_rectangular(t, t_start, duration, cell, E_field_magnitude):
    o = np.zeros_like(t)
    o += (delta_transmembrane_unit_step(t-t_start, cell) - delta_transmembrane_unit_step(t-t_start-duration, cell))
    o *= E_field_magnitude
    return o

def delta_transmembrane_rectangular_train(t, t_start, duration, count, period, cell, E_field_magnitude):
    o = np.zeros_like(t)
    for i in range(0, count):
        o += (delta_transmembrane_unit_step(t-t_start-(i*(period)), cell) - delta_transmembrane_unit_step(t-t_start-duration-(i*(period)), cell))
    o *= E_field_magnitude
    return o

def delta_transmembrane_rectangular_array(t, count, x, cell):
    o = np.zeros_like(t)
    # start, duration,
    for i in range(0, count):
        o += (delta_transmembrane_unit_step(t-x[i*3], cell) - delta_transmembrane_unit_step(t-x[(i*3)+1], cell)) * 10e6 * x[i*3 + 2]
    return o

def delta_transmembrane_trapezoid(t, t_start, t_rise, t_fall, duration, cell, E_field_magnitude):
    # https://lpsa.swarthmore.edu/LaplaceXform/FwdLaplace/LaplaceFuncs.html#Ramp
    # When composing a complex function from elementary functions, it is important to only use addition.
    # If you create a function by adding two functions, its Laplace Transform
    # is simply the sum of the Laplace Transform of the two function[s].
    # If you create a function by multiplying two functions in time,
    # there is no easy way to find the Laplace Transform of the resulting function.


    #ramp slope is 1: constructed from various ramps
    t_shifted = t - t_start
    out =  delta_transmembrane_unit_ramp(t_shifted, cell)*ustep(t_shifted)
    out -= delta_transmembrane_unit_ramp(t_shifted - t_rise, cell)*ustep(t_shifted - t_rise)
    out -= delta_transmembrane_unit_ramp(t_shifted - (duration - t_rise), cell)*ustep(t_shifted - (duration - t_rise))
    out += delta_transmembrane_unit_ramp(t_shifted, cell)*(t_shifted-duration)*ustep(t_shifted - (duration))
    # is the duration t1+rise perhaps?

    return out * E_field_magnitude
