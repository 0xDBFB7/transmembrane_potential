include("nagurka_fourier_lib.jl")

using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "../../")
tl = pyimport("transmembrane_lib")

py"""
import numpy as np
"""

virus = tl.default_virus(py"""np.array([])""")
host_cell = tl.default_host_cell(py"""np.array([])""")



using Gnuplot

pore_N0, pore_alpha, pore_q, pore_V_ep = tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep




# see https://discourse.julialang.org/t/enums-for-array-indexing/669/8
@enum svars begin
    iu0
    iu1

    iN_v
    iN_h

    ix0_v
    ix1_v

    ix0_h
    ix1_h
end
Base.to_index(s::svars) = Int(s+1)

# https://diffeqparamestim.sciml.ai/dev/methods/optimization_based_methods/
# might be useful

# or, more directly, https://diffeqparamestim.sciml.ai/dev/

# https://gcalderone.github.io/Gnuplot.jl/v1.1.0/

# http://julianlsolvers.github.io/Optim.jl/v0.9.3/algo/nelder_mead/

# not obvious to me how paramestim can have a free end time?

M = 3

t_f = 1e-6

function d_pore_density(V_m, N, N0, alpha, q, V_ep)
    k = (V_m/V_ep)^2

    return alpha * exp(k) * (1 - (N/N0)*exp(-q*k))
end;



function transmembrane_diffeq(d,s,O,t)
    """
    S is current state vectors.
    d is state vector derivtatives.
    O is a vector of parameters.
    """

    a = O[1:M]
    b = O[M+1:(2*M)]
    c = O[(2*M)+1:(2*M)+5]
    
    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = c

    X_t0 = 0.0 # override
    X_tf = 0.0
    d_X_t0 = 0.0

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    U = X_(t,p,a,b,m,t_f)
    d_U = d_X_(t,p,a,b,m,t_f)
    d_d_U = d_d_X_(t,p,a,b,m,t_f)

    T0 = t_f

    d[ ix0_v ] = s[ix1_v]
    d[ ix1_v ] = ((U0 / (T0*T0))*virus.alpha*d_d_U + (U0 / T0)*virus.beta*d_U + virus.gamma*U0*U - virus.phi*(X0 / T0)*s[ix1_v] - virus.xi*X0*s[ix0_v])/(X0 / (T0*T0))

    
    d[iN_v] = d_pore_density(x0_v, s[iN_v], pore_N0, pore_alpha, pore_q, pore_V_ep)

    d[iN_h] = d_pore_density(x0_h, s[iN_h], pore_N0, pore_alpha, pore_q, pore_V_ep)

end




#@gp t "with lines"
