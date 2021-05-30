from transmembrane_lib import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from gekko import GEKKO
import numpy as np
from scipy import signal
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

virus = Cell(np.float128(0.0001), np.float128(8000), np.float128(0.05), np.float128(300), np.float128(1e-8), np.float128(60), np.float128(20e-6), np.float128(14e-9), t1)

x0_v = m.Var(value=0)
x1_v = m.Var(value=0)
x2_v = m.Var(value=0)

t = m.Param(value=m.time)

u_step = np.zeros(nt,dtype=np.float128)
u_step[nt//2:] = 1
u_step[nt//2] = 0.5 #heaviside step
u0 = m.Param(value=u_step)
u1_ = np.diff(u_step, prepend=u_step[0])/((end/nt)/T0)
u1 = m.Param(value=u1_)
u2_ = np.diff(u1_, prepend=u1_[0])/((end/nt)/T0)
u2 = m.Param(value=u2_)

alpha_v = m.Const(virus.alpha)
beta_v = m.Const(virus.beta)
gamma_v = m.Const(virus.gamma)
phi_v = m.Const(virus.phi)
xi_v = m.Const(virus.xi)

SX0 = m.Const(X0)
SU0 = m.Const(U0)
ST0 = m.Const(T0)

m.Equation(x1_v==x0_v.dt())
m.Equation(x2_v==x1_v.dt())
m.Equation(x2_v == ((SU0 / (ST0**2))*alpha_v*u2 + (SU0 / ST0)*beta_v*u1 + gamma_v*SU0*u0 - phi_v*(SX0 / ST0)*x1_v - xi_v*SX0*x0_v)/(SX0 / (ST0**2)))

m.options.IMODE = 4 # dynamic simulation

m.solve(disp=True)

virus_output = np.array(x0_v.value)

plt.figure(1) # plot results
integration_method = virus_output*X0 / U0
# plt.plot(m.time*T0,,'r',marker='o')

plt.plot(m.time*T0,integration_method,'r',marker='o', label="GEKKO integrator")
# plt.plot(m.time*T0,x1_v.value,'g',marker='o')

convolution_method_output = convolve_output(u_step, virus, end/nt)
plt.plot(t1, convolution_method_output, label="lib_onvolution_method_output")

# shift(virus.step_response.astype('float32',casting='same_kind'), nt//2,cval=0.0)
plt.plot(t1, virus.step_response.astype('float32',casting='same_kind'), marker='o', label='lib_step_response')

system = signal.lti(np.float64([virus.a_3,virus.a_2,virus.a_1]),np.float64([virus.b_3,virus.b_2,virus.b_1]))

# _, yout, _ = signal.lsim2(system, U=np.float64(u_step), T=np.float64(t1)) # this isn't working for some reason
# plt.plot(t1, yout, label='lsim2')
# _, yout, _ = signal.lsim(system, U=np.float64(u_step), T=np.float64(t1))
# plt.plot(t1, yout, label='lsim')
# _, yout = signal.step(system, T=np.float64(t1))
# plt.plot(t1+end/2, yout, label='scipy step')

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
