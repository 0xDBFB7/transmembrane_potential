'''

N.B.

Contains the actual virus optimization routine.

'''



from nagurka_fourier_lib import *

import sys
sys.path.app('../')
from transmembrane_lib import *


virus = default_virus(np.array([]))
default_host_cell = default_host_cell(np.array([]))


M = 1



def X_v_():
    X_v = X_(t,P_BCs,a,b,M,t_f)
    d_X_v = d_X_(t,P_BCs,a,b,M,t_f)
    d_d_X_v = d_d_X_(t,P_BCs,a,b,M,t_f)


def cost_function(guess)
    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])


    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    t = np.linspace(epsilon, t_f, 60)


guess_initial = np.zeros(M*2 + 3)


# Tmin = minimize(cost_function, T, method="Powell", options={"disp":True}, callback=diagnostics, bounds=bounds, tol=1e-12).x
# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Powell", options={"disp":True}, bounds=bounds, callback=diagnostics,  tol=1e-12)
# Tmin = tubthumper(cost_function, T, stepsize=t_end/10, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]
