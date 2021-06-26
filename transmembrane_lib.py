from dataclasses import dataclass
from math import pi, sqrt, e, log, isclose, exp
# from scipy.optimize import curve_fit
try:
    import autograd.numpy as np  # Thinly-wrapped numpy
except:
    import numpy as np
from numpy import heaviside
# import h5py
from scipy.constants import epsilon_0, mu_0
from scipy.signal import convolve
from pytexit import py2tex

def ustep(v):

    # note that the definition of the discrete unit step can be the cause of differences
    # in outputs between solvers.
    # some use 0.5, some 1.
    return heaviside(v, 0.5)

# def safe_ustep(output, xes, h0=0.5):
#     # multiplying inf * 0 causes undef behav. This avoids that issue.
#     # also, more compatible with autograd.
#     return np.nan_to_num(output*heaviside(xes, 0.5))
#     # if(x > 0):
#     #     return output
#     # if(x == 0):
#     #     return h0
#     # if(x < 0):
#     #     return 0


@dataclass
class Cell:
    extracellular_conductivity: np.float64 # S/m
    extracellular_permittivity: np.float64 # relative
    intracellular_conductivity: np.float64 # S/m
    intracellular_permittivity: np.float64 # relative
    membrane_conductivity: np.float64 # S/m
    membrane_permittivity: np.float64 # relative
    cell_diameter: np.float64 # meters
    membrane_thickness: np.float64

    t: np.ndarray

    def __post_init__(self):
        e_o = self.extracellular_permittivity * epsilon_0 # S/m
        e_i = self.intracellular_permittivity * epsilon_0 #S/m
        e_m = self.membrane_permittivity * epsilon_0 #S/m
        R = self.cell_diameter / 2.0
        self.R = R

        l_o = self.extracellular_conductivity # S/m
        l_i = self.intracellular_conductivity #S/m
        l_m = self.membrane_conductivity #S/m

        d = self.membrane_thickness

        sub1 = (3.0 * (R**2.0) - 3.0 * d * R + d**2.0)
        sub2 = (3.0 * d * R - d**2.0)

        self.a_1 = 3.0 * d * l_o * ((l_i * (sub1)) + l_m*(sub2)) #eq.9a
        self.a_2 = 3.0 * d * ((l_i * e_o + l_o * e_i) * sub1 + (l_m * e_o + l_o * e_m) * sub2)
        self.a_3 = 3.0 * d * e_o * (e_i * (sub1) + e_m * sub2)

        self.b_1 = 2.0 * R**3.0 * (l_m + 2.0*l_o) * (l_m + 0.5 * l_i) + 2.0 * (R-d)**3.0 * (l_m - l_o) * (l_i - l_m)

        self.b_2 = 2.0 * R**3.0 * (l_i * (0.5 * e_m + e_o) + l_m * (0.5*e_i + 2.0*e_m + 2*e_o) + l_o * (e_i + 2.0 * e_m)) + (2.0 * (R - d)**3.0\
        * (l_i * (e_m - e_o) + l_m * (e_i - 2.0*e_m + e_o) - l_o * (e_i - e_m))) # is this truly a multiply, or a cross?


        self.b_3 = 2.0 * R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)

        # to check, use
        # py2tex("b_3 = 2.0 - R**3.0 * (e_m + 2.0*e_o) * (e_m + 0.5 * e_i) + 2.0 * (R-d)**3.0 * (e_m - e_o) * (e_i - e_m)")
        # https://quicklatex.com/

        # Kotnik variously use "step function" or the "unit step".
        self.tau_1 = tau_1_f(self.b_1, self.b_2, self.b_3)
        self.tau_2 = tau_2_f(self.b_1, self.b_2, self.b_3)


        self.step_response = delta_transmembrane_unit_step(self.t, self)

        # coefficients in front of the differential equation.
        # note the error in kotnik et al transfer function!
        self.alpha = ((self.R*self.a_3/self.b_3))
        self.beta = ((self.R*self.a_2/self.b_3))
        self.gamma = ((self.R*self.a_1/self.b_3))
        self.phi = ((self.b_2/self.b_3))
        self.xi = ((self.b_1/self.b_3))




def default_host_cell(t):
    # if(args):
    return Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)

def default_virus(t):
    #should be implemented but I haven't the heart to break things right now
    # if(args):
    #     t = args[0]
    # else:
    #     t = np.array(0)
    return Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8), np.float128(60), np.float128(50e-9), np.float128(14e-9),t)
    
'''

This is a verbatim implementation of

Kotnik T, Miklavčič D, Slivnik T.
Time course of transmembrane voltage induced by time-varying electric fields—a method for theoretical analysis and its application.
Bioelectrochemistry and Bioenergetics 1998;45:3–16.
https://doi.org/10.1016/S0302-4598(97)00093-7.

Kotnik et al also have number of other papers with minor variations that may be of some use

If a truly arbitrary
the equivalent of the discrete fourier transform for the Laplace transform appears to be the Z-transform.

Originally Kotnik's pre-made Laplace Transform-constructed waveforms were used.

-----------------------------------------

another setup is in

[1]Talele S, Gaynor P. Non-linear time domain model of electropermeabilization:
Response of a single cell to an arbitrary applied electric field.
Journal of Electrostatics 2007;65:775–84. https://doi.org/10.1016/j.elstat.2007.06.004.

At first glance it seems like the step response doesn't consider the internal permittivity:
(a usual three-layer structure isn't used)
but worry not, the third parameter is baked into /e'2

On the other hand, there's only one time constant here. Is that going to significantly affect things?

https://dsp.stackexchange.com/questions/17035
"So, apart from the first constant term,
the output signal y(t) is given by the convolution of the derivative
of the input signal with the system's step response."

"Again, apart from the constant first term, the output sequence is obtained by
convolving the first order difference of the input sequence with the system's step response."

"You can also find the impulse response from the step response by differentiation.
That means an alternative way to calculate the output is convolving with the derivative of the step response."

-----------------------------------------


Scipy lti has built-in functions for the convolution integral. Should replace this with that.

'''

def tau_1_f(b_1, b_2, b_3):
    #equation A6e in Kotnik
    return (2.0 * b_3) / (b_2 - np.sqrt(((b_2)**2.0) - ((4.0*b_1) * b_3)))

def tau_2_f(b_1, b_2, b_3):
    #equation A6f in Kotnik
    return (2.0 * b_3) / (b_2 + np.sqrt(((b_2)**2.0) - ((4.0*b_1) * b_3)))




def delta_transmembrane_unit_ramp(t, cell):
    # a9d, Kotnik 1998, unit ramp function response
    delta_phi_m_t = (cell.a_1 / cell.b_1) * t * ustep(t)
    phisub1 = (cell.a_2 / (2 * cell.b_1)) - ((cell.a_1*cell.b_2) / (2 * (cell.b_1 **2.0)))
    phisub2 = (((cell.a_1 * cell.b_3)/cell.b_1) + ((cell.a_2 * cell.b_2) / (2.0 * cell.b_1)) - ((cell.a_1 * (cell.b_2**2.0)) / (2.0 * (cell.b_1**2.0))) - cell.a_3)
    phisub3 = np.sqrt(b_2**2.0 - 4.0 * b_1 * b_3)

    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/cell.tau_1)) * ustep(t)
    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/cell.tau_2)) * ustep(t) # glitch?
    delta_phi_m_t *= cell.R
    # # unlike in symbolic math, the unit step doesn't actually prevent errors from occuring in the np.exp step.
    delta_phi_m_t = np.nan_to_num(delta_phi_m_t)

    # if(np.isscalar(t)):
    #     if(t < 0):
    #         delta_phi_m_t == 0
    # if(np.isscalar(t)):


    return delta_phi_m_t

def delta_transmembrane_unit_step(t, cell):

    # Kotnik variously use "step function" or the "unit step". I think these have different definitions?

    # a9d, Kotnik 1998, unit step function response
    delta_phi_m_t = (cell.a_3 / cell.b_3)
    phisub1 = (cell.a_1 / (2.0 * cell.b_1)) - (cell.a_3 / (2.0 * cell.b_3))
    phisub2 = ((cell.a_1 * cell.b_2) /  (2.0 * cell.b_1)) - cell.a_2 + ((cell.a_3 * cell.b_2) / (2.0 * cell.b_3))
    phisub3 = np.sqrt(cell.b_2**2.0 - 4.0 * cell.b_1 * cell.b_3)

    delta_phi_m_t += (phisub1 + (phisub2/phisub3)) * (1.0 - np.exp(-t/cell.tau_1))
    delta_phi_m_t += (phisub1 - (phisub2/phisub3)) * (1.0 - np.exp(-t/cell.tau_2))
    delta_phi_m_t *= cell.R

    # unlike in symbolic math, the unit step doesn't actually prevent errors from occuring in the np.exp step.
    # fixed by safe_ustep
    delta_phi_m_t = np.nan_to_num(delta_phi_m_t)
    # if(np.isscalar(t)):
    #     if(t < 0):
    #         delta_phi_m_t == 0

    # if(np.isscalar(t)):

    #beware: unit step replaced with (t > 0),

    return delta_phi_m_t * (t > 0)


def convolve_output(input, cell, dt):

    #again, either the step response or the input can be differentiated.
    input = np.diff(input, prepend=input[0])/dt

    #scipy signal automatically chooses an FFT-based method
    #not sure why I'm clipping to half, and I'm not entirely sure if that's valid.
    #cell.t
    return (convolve(cell.step_response,input) * dt)[0:len(cell.t)]


def total_waveform_energy(input,dt):
    return np.sum(np.square(input)*dt)


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


class gekko_sim_method:
    def __init__(self, virus, host_cell, control_input_u0, u1=None, u2=None, N=500, T0=1e-6, X0=1e-6, U0=0):
        # U0 scaling is not really needed
        self.m = GEKKO() # initialize gekko
        self.m.options.MAX_ITER = 200
        U0 = 1.0
        X0 = 1e-6
        end = tstop / T0
        self.m.time = np.linspace(0,end,N)

        # control_input_u0 =
        # gekko_input = np.float64(interpolant(np.float64(m.time)*T0))
        # u0 = m.Param(value=gekko_input)
        # u1_ = np.diff(gekko_input, prepend=gekko_input[0])/((end/nt))
        # u1 = m.Param(value=u1_)
        # u2_ = np.diff(u1_, prepend=u1_[0])/((end/nt)/T0)
        # u2 = m.Param(value=u2_)

        self.x0_v = m.Var(value=0)
        self.x1_v = m.Var(value=0)
        self.x2_v = m.Var(value=0)

        self.x0_h = m.Var(value=0)
        self.x1_h = m.Var(value=0)
        self.x2_h = m.Var(value=0)

        alpha_h = m.Const(host_cell.alpha)
        beta_h = m.Const(host_cell.beta)
        gamma_h = m.Const(host_cell.gamma)
        phi_h = m.Const(host_cell.phi)
        xi_h = m.Const(host_cell.xi)

        alpha_v = m.Const(virus.alpha)
        beta_v = m.Const(virus.beta)
        gamma_v = m.Const(virus.gamma)
        phi_v = m.Const(virus.phi)
        xi_v = m.Const(virus.xi)

        SX0 = m.Const(X0)
        SU0 = m.Const(U0)
        ST0 = m.Const(T0)

        m.Equation(self.x1_v==self.x0_v.dt())
        m.Equation(self.x2_v==self.x1_v.dt())
        m.Equation(self.x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*self.x1_v - xi_v*SX0*self.x0_v)/(SX0 / (ST0**2)))

        m.Equation(self.x1_h==self.x0_h.dt())
        m.Equation(self.x2_h==self.x1_h.dt())
        m.Equation(self.x2_h == ((SU0 / (ST0**2))*alpha_h*u2 + (SU0 / ST0)*beta_h*u1 + gamma_h*SU0*u0 - phi_h*(SX0 / ST0)*self.x1_h - xi_h*SX0*self.x0_h)/(SX0 / (ST0**2)))

        m.options.IMODE = 4 # dynamic simulation

        m.solve(disp=True) # solve
