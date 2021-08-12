'''

N.B.

Contains the actual virus optimization routine.

'''

from nagurka_membrane_extensions import *

virus = default_virus(np.array([]))
host_cell = default_host_cell(np.array([]))


M = 1


def get_output(guess):
    a = np.array(guess[0:M])
    b = np.array(guess[M:(2*M)])
    t_f = guess[2*M]


    X_t0 = 0.0
    X_tf = [2*M+1]
    d_X_t0 = 0.0
    d_X_tf = [2*M+2]
    d_d_X_t0 = [2*M+3]
    d_d_X_tf = [2*M+4]

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    t = np.linspace(epsilon, t_f, 300)

    virus = default_virus(t)
    host_cell = default_host_cell(t)


    U = X_(t,P_BCs,a,b,M,t_f)
    virus_output = U_to_X(U, t, virus)
    host_cell_output = U_to_X(U, t, host_cell)

    return U, virus_output, host_cell_output, t


def cost_function(guess):

    U, virus_output, host_cell_output, t = get_output(guess)

    v1 = np.sum(virus_output*virus_output)
    h1 = np.sum(host_cell_output*host_cell_output)
    u1 = np.sum(U*U)

    print(guess[0:M], guess[M:(2*M)], t[-1], v1, h1)

    return -v1 + h1


guess_initial = np.zeros(M*2 + 4)
guess_initial[2*M] = 0.0000001

bounds = [(-10, 10)]*(2*M) + [(0, 0.001)] + [(None, None)] + [(None, None)] + [(None, None)]
# Tmin = minimize(cost_function, guess_initial, method="Powell", options={"disp":True, "maxiter":1000}, bounds=bounds, tol=1e-14).x
tubthumper = basinhopping
minimizer_kwargs = dict(method="Powell", options={"disp":True, "maxiter":10}, bounds=bounds, tol=1e-14) #, tol=1e-12
Tmin = tubthumper(cost_function, guess_initial, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]

U, virus_output, host_cell_output, t = get_output(Tmin)
virus = default_virus(t)
host_cell = default_host_cell(t)

plt.plot(t, U)
plt.plot(t, virus_output*1e7)
plt.plot(t, host_cell_output*1e7)
plt.show()
