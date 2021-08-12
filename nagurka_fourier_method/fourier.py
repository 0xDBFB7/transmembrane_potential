'''

N.B.

Contains the actual virus optimization routine.

'''

from nagurka_membrane_extensions import *

virus = default_virus(np.array([]))
host_cell = default_host_cell(np.array([]))


M = 5


def get_output(guess):
    m = np.arange(1, M+1)
    a = np.array(guess[0:M], dtype=np.float128)
    b = np.array(guess[M:(2*M)], dtype=np.float128)
    t_f = 1e-6


    X_t0 = 0.0
    X_tf = guess[2*M+1]  / t_f
    d_X_t0 = guess[2*M+5]*(1/t_f)
    d_X_tf = guess[2*M+2]*(1/t_f)
    d_d_X_t0 = guess[2*M+3]*(1/(t_f**2.0))
    d_d_X_tf = guess[2*M+4]*(1/(t_f**2.0))

    # ic(X_tf)
    # ic(d_X_tf)
    # ic(d_d_X_t0)
    # ic(d_d_X_tf)


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    t = np.linspace(epsilon, t_f, 1000, dtype=np.float128)#300

    virus = default_virus(t)
    host_cell = default_host_cell(t)


    U = L_(t,a,b,M,t_f)
    virus_output = U_to_X(U, t, virus)
    host_cell_output = U_to_X(U, t, host_cell)

    return U, virus_output, host_cell_output, t


def cost_function(guess):

    U, virus_output, host_cell_output, t = get_output(guess)

    v1 = np.sum(virus_output[virus_output > 0])
    h1 = np.sum(host_cell_output[host_cell_output > 0])
    u1 = np.sum(U*U)

    # print(guess[0:M], guess[M:(2*M)])
    import os
    os.system('clear')
    print("t_f = ", guess[2*M])
    print("val = ", abs(h1/v1), v1, h1)

    # return -v1 + h1
    return abs(h1/v1) # settles invariably at 35.8 if t bound is low
    # return abs(abs(host_cell_output[-1])/abs(virus_output[-1]))


guess_initial = np.array((np.random.random(M*2 + 6, )*2 - 1.0), dtype=np.float128)
# guess_initial[0] =
guess_initial[2*M] = 10**(-np.random.random()*15)

bounds = [(-1000, 1000.0)]*(2*M) + [(1e-17, 1e-9)] + [(-10, 10)] + [(-100, 100)] + [(-100, 100)] + [(-100, 100)] + [(-10, 10)]

# Tmin = minimize(cost_function, guess_initial, method="Powell", options={"disp":True}, bounds=bounds).x #, "maxiter":1000
# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Nelder-Mead", options={"disp":True, "maxiter":100}, bounds=bounds) #, tol=1e-12
# Tmin = tubthumper(cost_function, guess_initial, minimizer_kwargs=minimizer_kwargs, disp=True)["x"]
#
# Tmin = guess_initial
# U, virus_output, host_cell_output, t = get_output(Tmin)
# virus = default_virus(t)
# host_cell = default_host_cell(t)
t_f = 1e-6
a = [0,0,0,0,1]
b = [0,0,0,0,0.25]
t = np.linspace(epsilon, t_f, 1000, dtype=np.float128)#300
U = t*0.0
m = np.arange(1, M+1)
for i in m:
    print(i)
    v = i*(2*pi / t_f)
    U += a[i-1]*np.cos(v*t)
# U = L_(t,a,b,M,t_f)
# plt.subplot(2, 1, 1)
plt.plot(t, U)
# plt.subplot(2, 1, 2)
# plt.plot(t, virus_output*1e7)
# plt.plot(t, host_cell_output*1e7)
plt.show()
