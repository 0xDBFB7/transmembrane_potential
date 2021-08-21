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

    iI_ep_v
    iI_ep_h
end
Base.to_index(s::svars) = Int(s)+1
Base.to_indices(I::Tuple) = (Int(I[1])+1, I[2])

# this darn redundancy needed because virus is global pystruct and slow mutable
struct cell_struct 
    alpha
    beta
    gamma
    phi
    xi
    membrane_permittivity
    membrane_thickness
    cell_diameter
end

mutable struct transmembrane_params
    # this makes it hard to use DiffEqParamEstim, but I can't see a way to pass t_f to diffeq otherwise.
    cell_v
    cell_h
    a
    b
    p
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
    # if(abs(V_m) > 2.0)
    #     V_m = 2.0
    # end
    k = (V_m/V_ep)^2
    return alpha * exp(k) * (1.0 - (N/N0)*exp(-q*k))
end;


function dVdt_tm_potential_pore_current(V_m, N, cell)
    """
    The term of the first derivative of the transmembrane potential
    caused by the current flow through the 
    """
    F = 96485.332 # units?
    T = 295.0
    R = 8.314

    r_m = 0.76e-9 # pore radius constant
    pore_solution_conductivity = 13.0 
    # normally 0.1 mS/cm to S/m, but this peaks out transmembrane at like 8 rather than 3
    # where does 1.3 S/m come from? that's pretty high...
    w0 = 2.65 # ?
    n = 0.15
    v_m = (V_m + epsilon) * (F/(R*T))

    i_ep_term_1 = (pi * (r_m^2) * pore_solution_conductivity * v_m * R * T) / (F * cell.membrane_thickness)

    i_ep_term_2_divisor_1 = exp(v_m)*((w0*exp(w0-n*v_m) - (n*v_m)) / (w0 - (n*v_m)))
    i_ep_term_2_divisor_2 = -((w0*exp(w0+n*v_m) + (n*v_m)) / (w0 + (n*v_m)))
    i_ep_term_2 = exp(v_m - 1) / (i_ep_term_2_divisor_1 + i_ep_term_2_divisor_2)

    i_ep = i_ep_term_1 * i_ep_term_2
    # I = C_m dv/dt
    # dv_dt = I/C_m
    # permittivity already has the eps0 in it 
    A = 4*pi*(cell.cell_diameter/2)^2
    C_m = cell.membrane_permittivity * A / cell.membrane_thickness 

    # C_m and A should maybe go in transmembrane_lib

    if(V_m >= 0.0)
        return -(i_ep * N / C_m)
    else
        return (i_ep * N / C_m)
    end

end

#=

When is the irreversible regime reached?

Retelj:

> "If a pore density of N = 10^14 m^−2 was reached at the pole of a membrane (the point where it
was the highest), the membrane was considered to be electroporated (threshold of significant, i.e., observable electroporation).
This value was taken from the model of DeBruin and Krassowska [49], which compared simulations with the experiments of Hibino et al. [1]."

Improved2002 suggest:

"Irreversible breakdown and cell rupture would, therefore, be the predicted result, for pores exceeding
the stability threshold radius r crit of 18 nm. In Fig. 1, both the peak energy and radius of the local maxima shift for a 0.2 V
transmembrane potential. The critical radius for stability reduces to about 5.8 nm."

> "First, as evident from Fig. 1, there is no barrier for 0.4 V. 
However, from experimental data, much higher membrane voltages of about 
1.0 V are required [28] for irreversible breakdown and membrane rupture."

DeBruin Modeling1999a use a constant pore density of 0.76 nm, and state that:

> "Including the effects of pore radius requires a substantial addition to the model that will be the subject of a future study."

=#


#=

Added an exp limiter to the transmembrane voltage,
representing the current flow through the pores.
the exp limiter could also be moved to the pore density term.

I have now moved this to the pore density term. 

~~This is going against my own views on standardization. 
I should use the validated pore current equations.~~
Done.

The problem with putting an * exp_limiter(iN_v)
on the first derivative of transmembrane potential 
is that the potential never decreases again.

The sharp rise of N causes a huge numerical instability.

An artificial value of conductivity has been used in the pore.

=#

function transmembrane_diffeq(d,s,params::transmembrane_params,t)
    """
    S is current state vectors.
    d is state vector derivtatives.
    O is a vector of parameters.
    """
    
    m = params.m
    M = params.M
    t_f = params.t_f

    @timeit to "us" begin
        s[iu0] = U = X_(t,params.p,params.a,params.b,m,t_f)
        s[iu1] = d_U = d_X_(t,params.p,params.a,params.b,m,t_f)
        s[iu2] = d_d_U = d_d_X_(t,params.p,params.a,params.b,m,t_f)
    end

    T0 = 1.0
    U0 = 1.0
    X0 = 1.0
    
    irreversible_threshold = 1e20
    irreversible_threshold_sharpness = 2.0 # sharper values seem to cause instabilities


    second_derivative_eq(cell, s1, s0) = (((U0 / (T0*T0))*cell.alpha*d_d_U + (U0 / T0)*cell.beta*d_U 
                                    + cell.gamma*U0*U - cell.phi*(X0 / T0)*s[s1] - cell.xi*X0*s[s0])/(X0 / (T0*T0)))

    exp_limiter(iN) = (exp(-(s[iN] / irreversible_threshold)^irreversible_threshold_sharpness))

    s[iI_ep_v] = (dVdt_tm_potential_pore_current(s[ix0_v], s[iN_v], params.cell_v) / T0) #? t0 right?
    s[iI_ep_h] = (dVdt_tm_potential_pore_current(s[ix0_h], s[iN_h], params.cell_h) / T0) #? t0 right?

    @timeit to "diffeq" begin
    d[ ix0_v ] = s[ix1_v] + s[iI_ep_v]
    d[ ix1_v ] =  second_derivative_eq(params.cell_v, ix1_v, ix0_v) 

    d[ ix0_h ] = s[ix1_h] + s[iI_ep_h]
    d[ ix1_h ] = second_derivative_eq(params.cell_h, ix1_h, ix0_h)
    end
    
    d[iN_v] = (d_pore_density(s[ix0_v], s[iN_v], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0) * exp_limiter(iN_v)
    d[iN_h] = (d_pore_density(s[ix0_h], s[iN_h], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0) * exp_limiter(iN_h)
    # @show s
    # @show d


    return
end




#@gp t "with lines"
