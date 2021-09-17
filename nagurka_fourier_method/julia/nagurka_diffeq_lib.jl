using OrdinaryDiffEq
using Gnuplot

include("nagurka_fourier_lib.jl")
include("pore_lib.jl")
include("cell_lib.jl")

using Printf
using TimerOutputs
to = TimerOutput()
disable_timer!(to)
# call reset_timer!() to reset

using ForwardDiff
using DoubleFloats

# try --color

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())


# julia --project=dev_nagurka/

# see https://discourse.julialang.org/t/enums-for-array-indexing/669/8
@enum svars begin
    iN_v
    iN_h

    ix0_v
    ix1_v

    ix0_h
    ix1_h

    ix2_v
    ix2_h

    iu0
    iu1
    iu2

    iI_ep_v
    iI_ep_h

    ilm_ep_v
    ilm_ep_h

    ialpha 
    ibeta 
    igamma 
    iphi
    ixi
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
    pore_solution_conductivity
    pore_model_enabled
    T0
    ufun
    d_ufun
    d_d_ufun
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


Test 


# irreversible_threshold = 1e40
# irreversible_threshold_sharpness = 2.0 # sharper values seem to cause instabilities
# exp_limiter(iN) = (exp(-(s[iN] / irreversible_threshold)^irreversible_threshold_sharpness))

=#


function main_x2_ode(d_d_U, d_U, U, s, cell, T0, l_m_ep, is1, is0)

    U0 = cell.xi / cell.gamma
    #important that this U0 doesn't change along with the conductivity!

    X0 = 1.0

    alpha, beta, gamma, phi, xi = electroporation_coefficients(cell, l_m_ep)

    # @printf("%f\n %f\n  %f\n  %f\n %f\n", ForwardDiff.value(alpha), ForwardDiff.value(beta), ForwardDiff.value(gamma), ForwardDiff.value(phi), ForwardDiff.value(xi))
    # @show gamma*U0*U (U0 / T0)*beta*d_U (U0 / (T0*T0))*alpha*d_d_U
    # alpha, beta, gamma, phi, xi = (cell.alpha, cell.beta, cell.gamma, cell.phi, cell.xi)
    return  (((U0 / (T0*T0))*alpha*d_d_U + (U0 / T0)*beta*d_U 
                                + gamma*U0*U - phi*(X0 / T0)*s[is1] - cell.xi*X0*s[is0])/(X0 / (T0*T0)))

end

function transmembrane_diffeq(d,s,params::transmembrane_params,t)
    """
    S is current state vectors.
    d is state vector derivtatives.
    O is a vector of parameters.
    """
    
    m = params.m
    M = params.M
    t_f = params.t_f
    T0 = params.T0

    s[iu0] = U = params.ufun(t)
    s[iu1] = d_U = params.d_ufun(t)
    s[iu2] = d_d_U = params.d_d_ufun(t)
    # @timeit to "us" begin
    #     s[iu0] = U = X_(t,params.p,params.a,params.b,m,t_f)
    #     s[iu1] = d_U = d_X_(t,params.p,params.a,params.b,m,t_f)
    #     s[iu2] = d_d_U = d_d_X_(t,params.p,params.a,params.b,m,t_f)
    # end

    # @timeit to "diffeq" begin
    
    # I_ep_v = electroporation_pore_current(s[ix0_v], s[iN_v], params.cell_v)
    # I_ep_v = 0.0
    # delta in conductivity due to the formed pores
    # if(abs(s[ix0_v]) > 0.0001)
    # l_m_ep_v = abs(I_ep_v / abs(s[ix0_v]))
    # else
    # end 
    # l_m_ep_v = 0.2 * (1 - exp(-s[iN_v] / 1e13))
    # s[ilm_ep_v] = l_m_ep_v
    # s[iI_ep_v] = I_ep_v
    
    # if(t < 0.0001*t_f)
    l_m_ep_v = 0.0
    # end

    # s[ ix2_v ] = main_x2_ode(d_d_U, d_U, U, s, params.cell_v, T0, l_m_ep_v, ix1_v, ix0_v)
    d[ ix1_v ] = s[ix2_v] 
    d[ ix0_v ] = s[ix1_v]
    
    # alpha, beta, gamma, phi, xi = electroporation_coefficients(params.cell_v, l_m_ep_v)
    # @show alpha, beta, gamma, phi, xi

    I_ep_h = electroporation_pore_current(s[ix0_h], s[iN_h], params.cell_h, params.pore_solution_conductivity)

    # if(abs(s[ix0_h]) > 0.0001)
    # l_m_ep_h = abs(I_ep_h / abs(s[ix0_h]))
    l_m_ep_h = current_to_conductivity(I_ep_h, s[ix0_h], params.cell_h.membrane_thickness, params.cell_h.cell_diameter / 2)
    # else
        
    # end
    
    # 

    if(abs(s[ix0_h]) < 0.01)
        l_m_ep_h = 0.0
    end

    if(l_m_ep_h > 1.3)
        l_m_ep_h = 1.3
    end

    s[ilm_ep_h] = l_m_ep_h
    s[iI_ep_h] = I_ep_h

    s[ ix2_h ] = main_x2_ode(d_d_U, d_U, U, s, params.cell_h, T0, l_m_ep_h, ix1_h, ix0_h)
    d[ ix1_h ] = s[ix2_h] 
    d[ ix0_h ] = s[ix1_h]


    alpha, beta, gamma, phi, xi = electroporation_coefficients(params.cell_h, l_m_ep_h)

    s[ialpha] = alpha
    s[ibeta] = beta
    s[igamma] = gamma
    s[iphi] = phi
    s[ixi] = xi
    # @show d
    # @show s

    # end
    # irreversible_threshold = 1e15
    # irreversible_threshold_sharpness = 2.0 # sharper values seem to cause instabilities
    # exp_limiter(iN) = (exp(-(s[iN] / irreversible_threshold)^irreversible_threshold_sharpness))

    # # this causes a lot of headache because of the large exponent - would be great to nondimensionalize this somehow, make it log perhaps
    # Adding a / T0 here makes the re-sealing time way unphysically fast - must not be correct.
    # not sure if it should be * T0
    # d[iN_v] = d_pore_density(s[ix0_v], s[iN_v], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep) * exp_limiter(iN_v)
    d[iN_h] = d_pore_density(s[ix0_h], s[iN_h], params.pore_N0, params.pore_alpha, params.pore_q, params.pore_V_ep)  #* exp_limiter(iN_h)

    return
end

function plot_solution(solution)

    formatstring = "with lines ls -1 dt 1 tit "
    @gp "set multiplot layout 5,2; set grid xtics ytics; set grid;"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) string(formatstring,"'u0'")
    @gp :- 2 solution.t getindex.(solution.u, Int(iu2)+1) string(formatstring,"'u2'")

    @gp :- 2 solution.t getindex.(solution.u, Int(iu1)+1) string(formatstring,"'u1'")
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_v)+1) string(formatstring,"'x0_v'")
    @gp :- 4 solution.t getindex.(solution.u, Int(ix0_h)+1) string(formatstring,"'x0_h'")
    @gp :- 5 solution.t (getindex.(solution.u, Int(iN_v)+1).- tl.pore_N0) string(formatstring,"'N_v'")
    @gp :- 6 solution.t (getindex.(solution.u, Int(iN_h)+1).- tl.pore_N0) string(formatstring,"'N_h'")
    @gp :- 7 solution.t (getindex.(solution.u, Int(ilm_ep_v)+1)) string(formatstring,"' \\lambda_{ep V}'")
    @gp :- 8 solution.t (getindex.(solution.u, Int(ilm_ep_h)+1)) string(formatstring,"' \\lambda_{ep H}'")
    @gp :- 9 solution.t (getindex.(solution.u, Int(iI_ep_v)+1)) string(formatstring,"' \$ I_{ep V}\$'")
    @gp :- 10 solution.t (getindex.(solution.u, Int(iI_ep_h)+1)) string(formatstring,"' \$ I_{ep H}\$'")


end

function schwan_steady_state_with_conductivities(E, G_m, cell)

    # eq. 7 from Dielectrophoresis and rotation of Cells, in the book 
    # G_m is the membrane conductance (S/m^2)
    # E is electric field intensity (V/m)
    return 1.5 * E * (cell.cell_diameter/2) / (1 + (cell.cell_diameter/2)*G_m*
                        ((1/cell.intracellular_conductivity) + 0.5*(1/cell.extracellular_conductivity))) 

end

#@gp t "with lines"
