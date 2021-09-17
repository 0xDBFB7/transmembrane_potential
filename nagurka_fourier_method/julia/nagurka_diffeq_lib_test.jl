include("nagurka_diffeq_lib.jl")

# BenchmarkTools.DEFAULT_PARAMETERS.samples = 30

# using Revise 
# this should probably be included before everything else

#using CUDA # just wrap initial_state_variables with cu()

#=

To do mem alloc analysis, 

julia --track-allocation=user
<run function> - precompiles
Profile.clear_malloc_data()
<run function>

https://github.com/JuliaApproximation/FastGaussQuadrature.jl apparently
is better if you can chose where to sample the thing, like an odesolver lets you.


Vsc shortcuts:
tree: ctl 0
editor: ctl 1
term: ctl k
close tab: ctl w
=#

# staticarrays.jl


@testset "Comparison" begin
    

    
    # @testset "time" begin
    #     d = zeros(length(instances(svars)))
    #     s = zeros(length(instances(svars)))
    #     t=epsilon
        
    #     single_timestep_time = @belapsed transmembrane_diffeq($d,$s,$params,t2) setup=(t2 = rand()) seconds=0.5
    #     @show single_timestep_time
    #     @show length(solution.t) * single_timestep_time
    # end
    
    # integrate(solution.t, getindex.(solution.u, Int(iu0)+1))
    # @btime integrate($solution.t, getindex.($solution.u, Int(iu0)+1))
    
end





function solve_response(params)

    tspan = (epsilon, params.t_f)
    initial_state_variables = (zeros(length(instances(svars)))).+epsilon #  convert(Array{BigFloat},zeros(length(instances(svars))))
    initial_state_variables[iN_v] = tl.pore_N0
    initial_state_variables[iN_h] = tl.pore_N0
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)
    solution = solve(prob, RadauIIA5(), dtmax = params.t_f / 200, maxiters= 1000000, dtmin=1e-20, 
                                                            progress = true, progress_steps = 100)
    #Tsit5
    return solution
end

function cell_step_test()

end


# @testset "Kotnik validation" begin
#     # Note: compare specify cell radius, not diameter, as 15e-6! 
#     cell_radius = 10e-6
#     compare_cell = tl.Cell((1.2), (80), (0.3), (80), (3e-7), (5), (cell_radius * 2), (5e-9), py"""np.array([])""")

#     cell_v = cell_h = py_cell_to_julia_struct(compare_cell)

#     end_time = 10e-6
#     T0 = 1e-6
#     # this is the "apparent" nondimensionalized end time to the algorithm
#     t_f = end_time / T0

#     edge_rise_time = 1e-9
#     k = 1 / (edge_rise_time / T0)
#     x0 = end_time / 100

#     E = 52000

#     peak = E * (cell_h.gamma / cell_h.xi)  # ((3/2) * (compare_cell.cell_diameter / 2)) # 52 kV/m * cell diameter scale factor

#     ufun(t) = logistic_curve( t, peak, k, x0)
#     d_ufun(t) = d_logistic_curve( t, peak, k, x0)
#     d_d_ufun(t) = d_d_logistic_curve( t, peak, k, x0)

    
#     pore_solution_conductivity = 0.6
#     params = transmembrane_params(cell_v, cell_h, nothing, nothing, nothing, t_f, nothing, nothing, tl.pore_N0, tl.pore_alpha, 
#                             tl.pore_q, tl.pore_V_ep, pore_solution_conductivity, false,  T0, ufun, d_ufun, d_d_ufun)
                            
#     solution = solve_response(params)

#     N_h_course = getindex.(solution.u, Int(iN_h)+1) # note: not subtracted from N0
    
#     plot_solution(solution)
    
#     compare_cell.compute_step_response(solution.t*T0)
    
#     G_m = 1.9 # membrane conductance

#     schwan_analytic = schwan_steady_state_with_conductivities(E, G_m, compare_cell)

#     @show schwan_analytic 

# end

@testset "Talele transmembrane potential validation" begin

    tal_ref_pore_density_interpolate, tal_ref_transmembrane_voltage_interpolate = import_talele_test_data()


    # Note: compare specify cell radius, not diameter, as 15e-6! 
    cell_radius = 15e-6
    compare_cell = tl.Cell((1.2), (80), (0.3), (80), (3e-7), (5), (cell_radius * 2), (5e-9), py"""np.array([])""")

    cell_v = cell_h = py_cell_to_julia_struct(compare_cell)

    end_time = 1e-6
    T0 = 1e-6
    # this is the "apparent" nondimensionalized end time to the algorithm
    t_f = end_time / T0

    edge_rise_time = 1e-9
    k = 1 / (edge_rise_time / T0)
    x0 = end_time / 100

    E = 52000

    peak = E * (cell_h.gamma / cell_h.xi)  # ((3/2) * (compare_cell.cell_diameter / 2)) # 52 kV/m * cell diameter scale factor

    ufun(t) = logistic_curve( t, peak, k, x0)
    d_ufun(t) = d_logistic_curve( t, peak, k, x0)
    d_d_ufun(t) = d_d_logistic_curve( t, peak, k, x0)

    
    pore_solution_conductivity = 0.6 
    params = transmembrane_params(cell_v, cell_h, nothing, nothing, nothing, t_f, nothing, nothing, tl.pore_N0, tl.pore_alpha, 
                            tl.pore_q, tl.pore_V_ep, pore_solution_conductivity, true, T0, ufun, d_ufun, d_d_ufun)
                            
    solution = solve_response(params)

    N_h_course = getindex.(solution.u, Int(iN_h)+1) # note: not subtracted from N0
    
    plot_solution(solution)

    @gp :- :GP2 1 solution.t transmembrane_voltage_interpolate(solution.t)
    @gp :- :GP2 2 solution.t pore_density_interpolate(solution.t)



    compare_cell.compute_step_response(solution.t*T0)
    
    G_m = 1.9 # membrane conductance

    schwan_analytic = schwan_steady_state_with_conductivities(E, G_m, compare_cell)

    @show schwan_analytic 

    # @gp :GP2 "set multiplot layout 2,1; set grid xtics ytics; set grid;"
    # @gp :- :GP2 1 solution.t getindex.(solution.u, Int(ix0_h)+1)
    # @gp :- :GP2 1 solution.t compare_cell.step_response * 52000

    # @gp :GP2 "set multiplot layout 5,1; set grid xtics ytics; set grid;"
    # @gp :- :GP2 1 solution.t getindex.(solution.u, Int(ialpha)+1)
    # @gp :- :GP2 2 solution.t getindex.(solution.u, Int(ibeta)+1)
    # @gp :- :GP2 3 solution.t getindex.(solution.u, Int(igamma)+1)
    # @gp :- :GP2 4 solution.t getindex.(solution.u, Int(iphi)+1)
    # @gp :- :GP2 5 solution.t getindex.(solution.u, Int(ixi)+1)
    
    @test_broken isapprox(getindex.(solution.u, Int(ix0_h)+1)[end], schwan_analytic)
    @test_broken isapprox(maximum(N_h_course), 3.334e13)

end

# @testset "DeBruin validation" begin
#     # DeBruin's formulation doesn't explicitly include the permittivity & dielectric relaxation of the cell, 
#     # only the membrane capacitance.
#     # Note: compare specify cell radius, not diameter, as 15e-6! 
#     cell_radius = 50e-6
#     E = 40000.0 # 400 V/cm to V/m
#     pore_solution_conductivity = 1.3
#     compare_cell = tl.Cell((5.0), (0.001), (0.455), (0.001), (3e-7), (5), (cell_radius * 2), (5e-9), py"""np.array([])""")

#     cell_v = cell_h = py_cell_to_julia_struct(compare_cell)

#     end_time = 10e-6
#     T0 = 1e-6
#     # this is the "apparent" nondimensionalized end time to the algorithm
#     t_f = end_time / T0

#     edge_rise_time = 1e-9
#     k = 1 / (edge_rise_time / T0)
#     x0 = end_time / 10

#     peak = E * (cell_h.gamma / cell_h.xi)  # ((3/2) * (compare_cell.cell_diameter / 2)) # 52 kV/m * cell diameter scale factor

#     ufun(t) = logistic_curve( t, peak, k, x0)
#     d_ufun(t) = d_logistic_curve( t, peak, k, x0)
#     d_d_ufun(t) = d_d_logistic_curve( t, peak, k, x0)

#     params = transmembrane_params(cell_v, cell_h, nothing, nothing, nothing, t_f, nothing, nothing, tl.pore_N0, tl.pore_alpha, 
#                             tl.pore_q, tl.pore_V_ep, pore_solution_conductivity, T0, ufun, d_ufun, d_d_ufun)
                            
#     solution = solve_response(params)

#     N_h_course = getindex.(solution.u, Int(iN_h)+1) # note: not subtracted from N0
    
#     plot_solution(solution)
# end