using DelimitedFiles
using Dierckx
using Test
using OrdinaryDiffEq
using Gnuplot


include("cell_lib.jl")
include("pore_lib.jl")
include("nagurka_fourier_lib.jl")


function pore_test_diffeq(d, s, p, t)
    d[1] = d_pore_density(p(t), s[1], tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep)
end

@testset "Talele pore formation rate comparison" begin
    """
    Since the 

    Uses the transmembrane voltages digitized from the plots in Talele et al 2007 fig. 2b.

    The peak of the pore is quite sensitive to the accuracy of the V_m digitization,
    so some deviation can be expected.
    """

    pore_density = readdlm("../../test_data/talele/talele_fig_2.csv", ',', skipstart=1, Float64)
    transmembrane_voltage = readdlm("../../test_data/talele/talele_fig_2_2.csv", ',', skipstart=1, Float64)
    
    microsecond = 1e-6
    pore_density_interpolate = Spline1D(pore_density[:,1] * microsecond, pore_density[:,2])
    transmembrane_voltage_interpolate = Spline1D(transmembrane_voltage[:,1] * microsecond, transmembrane_voltage[:,2])

    tspan = (epsilon, 1e-6)

    initial_state_variables = zeros(1)
    initial_state_variables[1] = tl.pore_N0
    prob = ODEProblem(pore_test_diffeq,initial_state_variables,tspan,transmembrane_voltage_interpolate)
    solution = solve(prob, Tsit5(), dtmax=1e-6/200, maxiters= 1000000, dtmin=1e-20, 
                                                            progress = true, progress_steps = 100)

    @gp :GP2 "set multiplot layout 2,1; set grid xtics ytics; set grid;"
    @gp :- :GP2 1 solution.t transmembrane_voltage_interpolate(solution.t)
    @gp :- :GP2 2 solution.t getindex.(solution.u, 1)
    @gp :- :GP2 2 solution.t pore_density_interpolate(solution.t)
    # @code_warntype    

    @test isapprox(pore_density_interpolate(1e-6), getindex.(solution.u, 1)[end], rtol=0.06)

end

# the problem with aborting the sim early based on N * pore area > membrane area
# is that it will mess up the cost function
# (because of the hard limit)
# https://diffeq.sciml.ai/stable/features/callback_functions/#Example-2:-Terminating-an-Integration