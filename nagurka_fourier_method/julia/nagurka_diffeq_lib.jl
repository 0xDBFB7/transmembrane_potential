using DifferentialEquations
using Gnuplot

include("nagurka_fourier_lib.jl")
include("pore_lib.jl")
include("cell_lib.jl")

using TimerOutputs
to = TimerOutput()
# disable_timer!(to)
# call reset_timer!() to reset






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
    
    irreversible_threshold = 1e22
    irreversible_threshold_sharpness = 2.0 # sharper values seem to cause instabilities


    second_derivative_eq(cell, s1, s0) = (((U0 / (T0*T0))*cell.alpha*d_d_U + (U0 / T0)*cell.beta*d_U 
                                    + cell.gamma*U0*U - cell.phi*(X0 / T0)*s[s1] - cell.xi*X0*s[s0])/(X0 / (T0*T0)))

    exp_limiter(iN) = (exp(-(s[iN] / irreversible_threshold)^irreversible_threshold_sharpness))

    # s[iI_ep_v] = d_d_V_ep(s[ix0_v], s[ix1_v], s[iN_v], params.cell_v)
    # s[iI_ep_h] = d_d_V_ep(s[ix0_h], s[ix1_h], s[iN_h], params.cell_h)

    

    @timeit to "diffeq" begin
    d[ ix0_v ] = s[ix1_v]
    d[ ix1_v ] = second_derivative_eq(params.cell_v, ix1_v, ix0_v) #+ s[iI_ep_v]

    d[ ix0_h ] = s[ix1_h]
    d[ ix1_h ] = second_derivative_eq(params.cell_h, ix1_h, ix0_h) #+ s[iI_ep_h]
    end
    
    # d[iN_v] = (d_pore_density(s[ix0_v], s[iN_v], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0) * exp_limiter(iN_v)
    # d[iN_h] = (d_pore_density(s[ix0_h], s[iN_h], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)/T0) * exp_limiter(iN_h)
    # @show s
    # @show d


    return
end




#@gp t "with lines"
