from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numpy import exp
from numpy import sqrt

'''
using the analytic soluion to u' u'' u = 0 provided by sympy_constraint,
test the output
'''

def analytic_control_solution(t, C1, C2):
    # without constraining the IC, c1 and c2 are different;
    # if u(t=0 ) = 0, C1=C2
    # C1 = - C2
    # base solution
    return (C1*np.exp(t*(-beta_h - np.sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))
    + C2*np.exp(t*(-beta_h + np.sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h)))
    # wait, why can this be imaginary?
    # what?
    # solution with ICS of u(0) = 0 u(1e-6)=1

    # return (np.exp(5.0e-7*(beta_h + 2.0*np.sqrt(-alpha_h*gamma_h + 0.25*beta_h**2))/alpha_h)*np.exp(t*(-beta_h + np.sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))/(np.exp(2.0e-6*np.sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h) - 1.0) +
    #     np.exp(5.0e-7*(beta_h + 2.0*np.sqrt(-alpha_h*gamma_h + 0.25*beta_h**2))/alpha_h)*np.exp(t*(-beta_h - np.sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))/(1.0 - np.exp(2.0e-6*np.sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h)))
    # just returns div 0
    # return (-2.0*alpha_h*exp(5.0e-7*(beta_h + 2.0*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2))/alpha_h)*exp(t*(-beta_h - sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))/(-beta_h*exp(2.0e-6*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h) + beta_h + 2.0*(-alpha_h*gamma_h + 0.25*beta_h**2)**0.5*exp(2.0e-6*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h) + 2.0*(-alpha_h*gamma_h + 0.25*beta_h**2)**0.5) + 2.0*alpha_h*exp(5.0e-7*(beta_h + 2.0*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2))/alpha_h)*exp(t*(-beta_h + sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))/(-beta_h*exp(2.0e-6*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h) + beta_h + 2.0*(-alpha_h*gamma_h + 0.25*beta_h**2)**0.5*exp(2.0e-6*sqrt(-alpha_h*gamma_h + 0.25*beta_h**2)/alpha_h) + 2.0*(-alpha_h*gamma_h + 0.25*beta_h**2)**0.5))
    #returns impulses at random intervals.
    # return 1*np.exp(t*(-beta_h + np.sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h)))


# a very strange and underwhelming result.
# apparently very short edges,
# very c
# at the largest scale, analytic_control_solution( with c1 and c2 -1 and 1 is a straight line
# with a slope that appears to be similar to what I've seen in previous simulations. However there's a fine structure
# of sub-femtosecond pulses.
# however, I don't understand why this convolution isn't agreeing with the expected value.
# we should get no response because u is 0  - oh it's because there's no second derivative value because it's between timesteps


t0 = 0
tstop = 1
t = np.linspace(t0, tstop, 1000000, dtype=np.float128)


host_cell = Cell(np.float128(0.3), np.float128(80), np.float128(0.3), np.float128(80), np.float128(1e-7), np.float128(5), np.float128(20e-6), np.float128(5e-9), t)
virus = Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-8), np.float128(60), np.float128(50e-9), np.float128(14e-9), t)

alpha_v = (virus.R*virus.a_1/virus.b_1)
beta_v = (virus.R*virus.a_2/virus.b_1)
gamma_v = (virus.R*virus.a_3/virus.b_1)
phi_v = (virus.b_2/virus.b_1)
xi_v = (virus.b_3/virus.b_1)


alpha_h = (host_cell.R*host_cell.a_1/host_cell.b_1)
beta_h = (host_cell.R*host_cell.a_2/host_cell.b_1)
gamma_h = (host_cell.R*host_cell.a_3/host_cell.b_1)
phi_h = (host_cell.b_2/host_cell.b_1)
xi_h = (host_cell.b_3/host_cell.b_1)


print(virus.R*virus.a_1)
print(host_cell.R*host_cell.a_2)
print(host_cell.R*host_cell.a_3)

 # differential equation is wrong - must be just after gekko_.
print(beta_h)
print(-4*alpha_h*gamma_h + beta_h**2)

control_input = analytic_control_solution(t, 1, 1)
control_input /= np.max(control_input)


plt.figure(1)
plt.plot(t, control_input,marker='o')
plt.show()
plt.figure(2)
# plt.plot(t, convolve_output(control_input, host_cell, dt))
# plt.plot(t, convolve_output(control_input, virus, dt))

print((np.max(convolve_output(control_input, host_cell, dt)) - np.min(convolve_output(control_input, host_cell, dt)))
            /(np.max(convolve_output(control_input, virus, dt)) - np.min(convolve_output(control_input, virus, dt))))

# plt.show()
