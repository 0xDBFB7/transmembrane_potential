from transmembrane_lib import *
from scipy.integrate import odeint

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 20e-6, 5e-9, np.array([]))
virus = Cell(0.3, 80, 0.005, 30, 1e-8, 60, 50e-9, 14e-9, np.array([]))

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


def analytic_control_solution(C1, C2):
    # without constraining the IC, c1 and c2 are different;
    # if u(t=0 ) = 0, C1=C2
    # C1 = - C2
     C1*exp(t*(-beta_h - sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h))
    + C2*exp(t*(-beta_h + sqrt(-4*alpha_h*gamma_h + beta_h**2))/(2*alpha_h)))
