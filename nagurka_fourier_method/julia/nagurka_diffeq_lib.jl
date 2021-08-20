using DifferentialEquations
using Gnuplot
include("nagurka_fourier_lib.jl")

using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "../../")
tl = pyimport("transmembrane_lib")

py"""
import numpy as np
"""

using TimerOutputs
to = TimerOutput()
disable_timer!(to)
# call reset_timer!() to reset


virus = tl.default_virus(py"""np.array([])""")
host_cell = tl.default_host_cell(py"""np.array([])""")





# see https://discourse.julialang.org/t/enums-for-array-indexing/669/8
@enum svars begin
    iN_v
    iN_h

    ix0_v
    ix1_v

    ix0_h
    ix1_h

    iu0
    iu1
    iu2
end
Base.to_index(s::svars) = Int(s)+1
Base.to_indices(I::Tuple) = (Int(I[1])+1, I[2])

struct cell_struct # this darn redundancy needed because virus is global pystruct and slow mutable
    alpha
    beta
    gamma
    phi
    xi
end

mutable struct transmembrane_params
    # this makes it hard to use DiffEqParamEstim, but I can't see a way to pass t_f to diffeq otherwise.
    cell_v
    cell_h
    O
    t_f 
    M 
    m
    pore_N0
    pore_alpha
    pore_q
    pore_V_ep
end


# https://diffeqparamestim.sciml.ai/dev/methods/optimization_based_methods/
# might be useful

# or, more directly, https://diffeqparamestim.sciml.ai/dev/

# https://gcalderone.github.io/Gnuplot.jl/v1.1.0/

# http://julianlsolvers.github.io/Optim.jl/v0.9.3/algo/nelder_mead/

# not obvious to me how paramestim can have a free end time?


function d_pore_density(V_m, N, N0, alpha, q, V_ep)
    k = (V_m/V_ep)^2

    return alpha * exp(k) * (1.0 - (N/N0)*exp(-q*k))
end;

function transmembrane_diffeq(d,s,params::transmembrane_params,t)
    """
    S is current state vectors.
    d is state vector derivtatives.
    O is a vector of parameters.
    """
    
    @timeit to "views" begin
    m = params.m
    M = params.M
    t_f = params.t_f


    @timeit to "us" begin
        s[iu0] = U = X_(t,p,a,b,m,t_f)
        s[iu1] = d_U = d_X_(t,p,a,b,m,t_f)
        s[iu2] = d_d_U = d_d_X_(t,p,a,b,m,t_f)
    end

    T0 = 1.0
    U0 = 1.0
    X0 = 1.0

    second_derivative_eq(cell, s1, s0) = ((U0 / (T0*T0))*cell.alpha*d_d_U + (U0 / T0)*cell.beta*d_U 
                                    + cell.gamma*U0*U - cell.phi*(X0 / T0)*s[s1] - cell.xi*X0*s[s0])/(X0 / (T0*T0))

    @timeit to "diffeq" begin
    d[ ix0_v ] = s[ix1_v]
    d[ ix1_v ] = second_derivative_eq(params.cell_v, ix1_v, ix0_v)
    d[ ix0_h ] = s[ix1_h]
    d[ ix1_h ] = second_derivative_eq(params.cell_h, ix1_h, ix0_h)
    end


    d[iN_v] = d_pore_density(s[ix0_v], s[iN_v], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0

    d[iN_h] = d_pore_density(s[ix0_h], s[iN_h], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0

    

end




#@gp t "with lines"
