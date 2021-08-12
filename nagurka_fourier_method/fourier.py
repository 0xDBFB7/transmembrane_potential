'''

N.B.

Contains the actual virus optimization routine.

'''

from nagurka_membrane_extensions import *

"""
Both bottom out at about 2.7, no obvious difference.
"""
def new_virus(t):
    return Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-7), np.float128(5), np.float128(50e-9), np.float128(14e-9),t)




M = 15

def get_output(guess):
    m = np.arange(1, M+1)
    a = np.array(guess[0:M], dtype=np.float128) #* m**2.0
    b = np.array(guess[M:(2*M)], dtype=np.float128)
    t_f = 1e-7


    X_t0 = 0.0
    X_tf = guess[2*M+1]
    d_X_t0 = 0.0 #guess[2*M+5]*(1/t_f)
    d_X_tf = guess[2*M+2]*(1/t_f)
    d_d_X_t0 = guess[2*M+3]*(1/(t_f**2.0))
    d_d_X_tf = guess[2*M+4]*(1/(t_f**2.0))

    # ic(X_tf)
    # ic(d_X_tf)
    # ic(d_d_X_t0)
    # ic(d_d_X_tf)


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    t = np.linspace(epsilon, t_f, 5000, dtype=np.float128)#300 before - really running into resolution issues

    virus = default_virus(t)
    host_cell = default_host_cell(t)


    U = L_(t,a,b,M,t_f)
    virus_output = U_to_X(U, t, virus)
    host_cell_output = U_to_X(U, t, host_cell)

    return U, virus_output, host_cell_output, t


def cost_function(guess):

    U, virus_output, host_cell_output, t = get_output(guess)

    # v1 = np.sum(virus_output[virus_output > 0])
    # h1 = np.sum(host_cell_output[host_cell_output > 0])
    # v1 = np.sum(virus_output*virus_output)
    # h1 = np.sum(host_cell_output*host_cell_output)

    v1 = np.abs(np.sum(virus_output))
    h1 = np.abs(np.sum(host_cell_output))


    u1 = np.sum(U*U)

    # print(guess[0:M], guess[M:(2*M)])
    import os
    os.system('clear')
    print("t_f = ", guess[2*M])
    print("val = ", abs(h1/v1), v1, h1)

    return -v1 + h1
    # return abs(h1/v1) # settles invariably at 35.8 if t bound is low
    # return abs(abs(host_cell_output[-1])/abs(virus_output[-1]))


guess_initial = np.array((np.random.random(M*2 + 5, )*2 - 1.0), dtype=np.float128)
# guess_initial[0] =
guess_initial[2*M] = 10**(-np.random.random()*15)

bounds = [(-1000, 1000.0)]*(2*M) + [(1e-10, 1e-4)] + [(-10, 10)] + [(-100, 100)] + [(-100, 100)] + [(-100, 100)] #+ [(-10, 10)]

Tmin = minimize(cost_function, guess_initial, method="Nelder-Mead", options={"disp":True, "maxiter":10000}, bounds=bounds).x #, "maxiter":1000
# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Nelder-Mead", options={"disp":True, "maxiter":100}, bounds=bounds) #, tol=1e-12
# Tmin = tubthumper(cost_function, guess_initial, minimizer_kwargs=minimizer_kwargs, disp=True)["x"]

U, virus_output, host_cell_output, t = get_output(Tmin)
U = U/np.max(U)
virus = default_virus(t)
host_cell = default_host_cell(t)
plt.subplot(2, 1, 1)
plt.plot(t, U)
plt.subplot(2, 1, 2)
plt.plot(t, virus_output*1e7)
plt.plot(t, host_cell_output*1e7)
v1 = np.sum(np.abs(virus_output))
h1 = np.sum(np.abs(host_cell_output))
plt.savefig(f"plots/{ abs(h1/v1)}.png")
plt.show()
