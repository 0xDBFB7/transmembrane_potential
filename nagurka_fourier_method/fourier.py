'''

N.B.

Contains the actual virus optimization routine.

'''
import dill
import pickle
import os
from nagurka_membrane_extensions import *

"""
Both bottom out at about 2.7, no obvious difference.
"""
def new_virus(t):
    return Cell(np.float128(0.3), np.float128(80), np.float128(0.005), np.float128(30), np.float128(1e-7), np.float128(5), np.float128(50e-9), np.float128(14e-9),t)

"""
A good intro to control parameterization is in
A Rigorous Global Optimization Algorithm for Problems with Ordinary Differential Equations

it seems like Lagrange polynomials are usually used for this sort of thing

"""

"""
https://en.wikipedia.org/wiki/Adaptive_coordinate_descent might be interesting
"""

M = 30 # number of fourier terms

#input_amplitude = 1e8#

vir_w = 1
"""
A watched nelder-mead never converges

depending on the time scale involved,
Has to be scaled by 1/20 initially to push the N into the nonlinear - otherwise the
host N count is beyond float precision.

Initially, transmembrane voltage clipping was used to prevent N from running off the edge.
(in reality, of course, the current flow due to the high number of pores would limit the voltage).
However, a much more effective technique was to scale the
"""

"""
There seems to be a weird linear increase in N that I don't think should be there.
"""

def get_output(guess):
    m = np.arange(1, M+1)
    a = np.array(guess[0:M], dtype=np.float128) #* m**2.0
    b = np.array(guess[M:(2*M)], dtype=np.float128)
    t_f = abs(guess[2*M])

    # input_amplitude = abs(guess[2*M]) * 1e7
    input_amplitude = 1


    # next, try adding a 1e7 scaling to a and b
    # and lengthening the time scale even more
    # maybe fiddle with the N integration algorithm? (is odeint really necessary?)

    ts = int(((2*pi*M))*7) # number of time steps

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

    t = np.linspace(epsilon, t_f, ts, dtype=np.float128)#300 before - really running into resolution issues

    virus = new_virus(t)
    host_cell = default_host_cell(t)

    U = X_(t,P_BCs,a,b,M,t_f)
    print(np.max(U))
    U /= np.max(np.abs(U))
    virus_output = U_to_X(U, t, virus) * input_amplitude
    host_cell_output = U_to_X(U, t, host_cell) * input_amplitude

    w = np.max(np.abs(host_cell_output))/3 + np.max(np.abs(virus_output))/3
    virus_output /= w/vir_w
    host_cell_output /= w
    # assumes the parameters are identical for both membranes
    Nsq_virus = integrate_pore_density(t, virus_output, pore_N0, pore_alpha, pore_q, pore_V_ep)
    Nsq_host_cell = integrate_pore_density(t, host_cell_output, pore_N0, pore_alpha, pore_q, pore_V_ep)

    return U, virus_output, host_cell_output, Nsq_virus, Nsq_host_cell, t, ts


def cost_function(guess):

    U, virus_output, host_cell_output, Nsq_virus, Nsq_host_cell, t, ts = get_output(guess)
    # v1 = np.sum(virus_output[virus_output > 0])
    # h1 = np.sum(host_cell_output[host_cell_output > 0])
    # v1 = np.sum(virus_output*virus_output)
    # h1 = np.sum(host_cell_output*host_cell_output)
    #
    # v1 = np.abs(np.sum(virus_output))
    # h1 = np.abs(np.sum(host_cell_output))

    #0.25 threshold almost works but of course it can jump quickly
    # honestly this works way better than you might expect
    #previously had np.where after sum
    # v1 = np.sum(([virus_output*1e7 > 0.5*np.max(virus_output*1e7)])) + epsilon
    # h1 = np.sum((np.logical_or([host_cell_output*1e7 > 0.5*np.max(host_cell_output*1e7)],
    #                     [host_cell_output*1e7 < 0.5*np.min(host_cell_output*1e7)]))) + epsilon

    # v1 = np.sum(([virus_output*1e7 > 0.1])) + epsilon  + np.sum(np.abs(virus_output))*ts
    # h1 = np.sum((np.logical_or([host_cell_output*1e7 > 0.1],
    #                     [host_cell_output*1e7 < -0.1]))) + epsilon + np.sum(np.abs(host_cell_output))*ts



    #0.25 threshold almost works but of course it can jump quickly # M=40
    # np. where isn't actually right, not what I was going for
    # v1 = np.sum(np.where([virus_output*1e7 > host_cell_output*1e7])) + epsilon
    # h1 = np.sum(np.where([host_cell_output*1e7 > virus_output*1e7]
    #                     )) + epsilon

    #can also reverse order of sum abs
    # v1 = np.sum(np.abs(virus_output*1e7 - 1.0)) # note with these the val is reversed (needs to increase)
    # h1 = np.sum(np.abs(host_cell_output*1e7 - 0.0))
    # generated 212.75299946576507.png

    # v1 = np.sum(([np.abs(virus_output*1e7) > 0.25]))
    # h1 = np.sum(([np.abs(host_cell_output*1e7) > 0.25]))

    # u1 = np.sum(U*U) #needs

    # v1 = np.sum(np.abs(virus_output*1e7 - 1.0) / (ts))
    # h1 = np.sum(np.abs(host_cell_output*1e7)/(ts))
    #
    # v2 = np.sum(np.abs(np.diff(virus_output))/(np.max(virus_output)-np.min(virus_output)))*0.01
    # h2 = np.sum(np.abs(np.diff(host_cell_output))/(np.max(host_cell_output)-np.min(host_cell_output)))*0.01

    # using max might be a bit unstable. maybe sum is better?
    # not really what we want
    v1 = np.sum(Nsq_virus)#[-1]
    h1 = np.sum(Nsq_host_cell)#[-1]


    # print(guess[0:M], guess[M:(2*M)])
    os.system('clear')
    print("t_f = ", abs(guess[2*M]))
    print("val = ", abs(h1/v1), v1, h1)#, v2, h2)


    return -v1 + h1 #+ v2 + h2 #-v1
    # return abs(h1/v1) # settles invariably at 35.8 if t bound is low
    # return abs(abs(host_cell_output[-1])/abs(virus_output[-1]))



# guess_initial = np.array((np.random.random(M*2 + 5, )*2 - 1.0), dtype=np.float128)
guess_initial = np.ones(M*2 + 5, dtype=np.float128)

try:
    with open('data.pickle', 'rb') as f:
        guess_initial = pickle.load(f)
except:
    pass

# vir_w = 1.0
# guess_initial[2*M] = 10**(-np.random.random()*15) #time
# guess_initial[2*M] = 1.0

# bounds = [(-1000, 1000.0)]*(2*M) + [(1e-12, 1e-4)] + [(-10, 10)] + [(-100, 100)] + [(-100, 100)] + [(-100, 100)] #+ [(-10, 10)]



# filename = 'globalsave.pkl'
# try:
#     dill.load_session(filename)
#     print()
#     print("#"*10 + "LOADED PREVIOUS SESSION" + "#"*10)
#     print()
# except:
#     Tmin = minimize(cost_function, guess_initial, method="Nelder-Mead", options={"disp":True, "maxiter":500}).x #, "maxiter":1000
#
#     dill.dump_session(filename)

Tmin = minimize(cost_function, guess_initial, method="Nelder-Mead", options={"disp":True, "maxiter":10000}).x #, "maxiter":1000


# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Nelder-Mead", options={"disp":True, "maxiter":10}) #, bounds=bounds, tol=1e-12
# Tmin = tubthumper(cost_function, guess_initial, minimizer_kwargs=minimizer_kwargs, disp=True)["x"]


with open('data.pickle', 'wb') as f:
    pickle.dump(Tmin, f)

U, virus_output, host_cell_output, Nsq_virus, Nsq_host_cell, t, ts = get_output(Tmin)

virus = new_virus(t)
host_cell = default_host_cell(t)
plt.subplot(3, 1, 1)
plt.plot(t, U)
plt.subplot(3, 1, 2)
plt.plot(t, virus_output)
plt.plot(t, host_cell_output)
v1 = np.abs(np.sum(virus_output))
h1 = np.abs(np.sum(host_cell_output))
plt.subplot(3, 1, 3)
plt.plot(t, Nsq_virus)
plt.plot(t, Nsq_host_cell)
plt.savefig(f"plots/{ abs(h1/v1)}.png")
plt.show()




# plt.plot(t, np.cumsum((([virus_output*1e8 > 0.25])))) # old way to integrate pore count
# plt.plot(t, np.cumsum((([host_cell_output*1e8 > 0.25]))))
