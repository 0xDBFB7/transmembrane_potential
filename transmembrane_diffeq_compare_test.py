from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO
import numpy as np

from scipy.ndimage.interpolation import shift


m = GEKKO() # initialize gekko
m.options.MAX_ITER = 30
nt = 1000
end = 1e-1

T0 = 1.0
U0 = 1.0
X0 = 1.0
m.time = np.linspace(0, end/T0, nt, dtype=np.float128)
t1 = np.linspace(0, end, nt, dtype=np.float128)

virus = Cell(np.float128(0.0001), np.float128(80), np.float128(0.05), np.float128(30), np.float128(1e-8), np.float128(60000), np.float128(50e-9), np.float128(14e-9), t1)

#override to check discretization error
virus.R = np.float128(1.0)
# virus.a_1 = virus.b_1 = np.float128(1.0)
# virus.a_2 = virus.b_2 = np.float128(1)
# virus.a_3 = virus.b_3 = np.float128(1e-4)
# virus.tau_1 = tau_1_f(virus.b_1, virus.b_2, virus.b_3)
# virus.tau_2 = tau_2_f(virus.b_1, virus.b_2, virus.b_3)

m.options.RTOL = 1e-1
m.options.OTOL = 1e-1

virus.step_response = delta_transmembrane_unit_step(virus.t, virus)

 # differential equation is wrong - must be just after gekko_.

control_input = np.sin(t1/((end/30.0)))



# Variables
x0_v = m.Var(value=0)
x1_v = m.Var()
# x2_v = m.Var()

t = m.Param(value=m.time)

# u0 = m.Var()
# u1 = m.Var()
# m.Equation(u1==u0.dt())
# u2 = m.Var()
# m.Equation(u2==u1.dt())

# m.Equation(u0 == 1)
# the gaussian pulse below works much better.
# m.Equation(u0 == m.sin(t/((end/30)/T0))) # for simulation
# m.Equation(u0 == m.exp(-((((t*T0)-((end*T0))/2.0))**2.0)/(2.0*(((((end*T0)/10.0))**2.0))))) # for simulation

u_step = np.zeros(nt,dtype=np.float128)
u_step[nt//2:] = 1
u_step[nt//2] = 0.5 #heaviside step
u0 = m.Param(value=u_step)
u1_ = np.diff(u_step, prepend=u_step[0])/((end/nt)/T0)
u1 = m.Param(value=u1_)
u2_ = np.diff(u1_, prepend=u1_[0])/((end/nt)/T0)
u2 = m.Param(value=u2_)


alpha_v = m.Const((U0 / (T0**2))*(virus.R*virus.a_3/virus.b_3))
beta_v = m.Const((U0 / T0)*(virus.R*virus.a_2/virus.b_3))
gamma_v = m.Const((virus.R*virus.a_1/virus.b_3))

phi_v = m.Const((X0 / T0)*(virus.b_2/virus.b_3))
xi_v = m.Const(X0*(virus.b_1/virus.b_3))





# phi_v = m.Const(0)
# xi_v = m.Const(0)

# alpha_v = m.Const(virus.R*virus.a_1/virus.b_1)
# beta_v = m.Const(virus.R*virus.a_2/virus.b_1)


SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
# m.Equation(x2_v==x1_v.dt())

print(alpha_v.value)
print(beta_v.value)
print(gamma_v.value)
print(phi_v.value)
print(xi_v.value)


#check coefficients
#should have the form of a  charging capacitor

m.Equation(x1_v.dt()==(alpha_v*u2 + beta_v*u1 + gamma_v*SU0*u0 - phi_v*x0_v.dt() - xi_v*x0_v))

m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True)

virus_output = np.array(x0_v.value)

plt.figure(1) # plot results
integration_method = virus_output*X0 / U0
# plt.plot(m.time*T0,,'r',marker='o')

plt.plot(m.time*T0,integration_method,'r',marker='o')
# plt.plot(m.time*T0,x1_v.value,'g',marker='o')

# convolution_method_output = convolve_output(control_input, virus, end/nt)
# plt.plot(t1, convolution_method_output)



plt.plot(t1, shift(virus.step_response.astype('float32',casting='same_kind'), nt//2,cval=0.0),marker='o', label='lib_step_response')



from scipy import signal
system = signal.lti(np.float64([virus.a_3,virus.a_2,virus.a_1]),np.float64([virus.b_3,virus.b_2,virus.b_1]))


# signal.tf2ss([virus.a_1,virus.a_2,virus.a_3],[virus.b_1,virus.b_2,virus.b_3])

_, yout, _ = signal.lsim2(system, U=np.float64(u_step), T=np.float64(t1))
plt.plot(t1, yout, label='lsim2')
_, yout, _ = signal.lsim(system, U=np.float64(u_step), T=np.float64(t1))
plt.plot(t1, yout, label='lsim')
_, yout = signal.step(system, T=np.float64(t1))
plt.plot(t1+end/2, yout, label='scipy step')

plt.legend()

plt.show()

#
#
#
# from scipy.integrate import quad
# def integrand(x, a, b):
#     return a*x**2 + b
#
# a = 2
# b = 1
# I = quad(integrand, 0, 1, args=(a,b))
