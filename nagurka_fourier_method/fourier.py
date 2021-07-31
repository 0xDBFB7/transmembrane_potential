from naguka_fourier_lib import *






t_f = 1e-6


a = np.array([0.0])
b = np.array([1.0])

M = a.shape[0]

t = np.linspace(0, t_f, 100)

virus = default_virus(t)
host_cell = default_host_cell(t)

X_t0 = 1.0
X_tf = 2.0
d_X_t0 = 0.0
d_X_tf = 0.0
d_d_X_t0 = 0.0
d_d_X_tf = 0.0



p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
output_course = P_(t,p,M) + L_(t,a,b,M,t_f)

plt.plot(t, output_course)
plt.show()


# Tmin = minimize(cost_function, T, method="Powell", options={"disp":True}, callback=diagnostics, bounds=bounds, tol=1e-12).x
# tubthumper = basinhopping
# minimizer_kwargs = dict(method="Powell", options={"disp":True}, bounds=bounds, callback=diagnostics,  tol=1e-12)
# Tmin = tubthumper(cost_function, T, stepsize=t_end/10, minimizer_kwargs=minimizer_kwargs, disp=True, niter_success=2)["x"]
