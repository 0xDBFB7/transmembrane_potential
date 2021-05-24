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

    return -C2*exp(t*(-host_cell.beta - sqrt(-4*host_cell.alpha*host_cell.gamma + host_cell.beta**2))/(2*host_cell.alpha)) + C2*exp(t*(-host_cell.beta + sqrt(-4*host_cell.alpha*host_cell.gamma + host_cell.beta**2))/(2*host_cell.alpha))

    # solution with ICS of u(0) = 0 u(1e-6)=1



t0 = 0
tstop = 1e-8
n = 10000
dt = tstop / n
t = np.linspace(t0, tstop, n, dtype=np.float128)


host_cell = default_host_cell(t)
virus = default_virus(t)


control_input = analytic_control_solution(t, -1, 1)

control_input /= np.max(control_input)

first_derivative = (np.diff(control_input)/dt)
second_derivative = np.diff(first_derivative)/dt

t_ = t[:-2]
control_input = control_input[:-2]
first_derivative = first_derivative[:-1]

sum = host_cell.alpha * second_derivative + host_cell.beta * first_derivative + host_cell.gamma * control_input

# plt.figure(1)
# plt.plot(t, host_cell.gamma * control_input,marker='o')
# plt.figure(2)
#
# plt.plot(host_cell.beta * first_derivative ,marker='o')
# plt.figure(3)
# plt.plot(host_cell.alpha * second_derivative ,marker='o')
# plt.figure(4)
plt.plot(t, convolve_output(control_input, host_cell, dt))
# plt.plot(t_, sum,marker='o')

plt.show()



print((np.max(convolve_output(control_input, host_cell, dt)) - np.min(convolve_output(control_input, host_cell, dt)))
            /(np.max(convolve_output(control_input, virus, dt)) - np.min(convolve_output(control_input, virus, dt))))

# plt.show()
