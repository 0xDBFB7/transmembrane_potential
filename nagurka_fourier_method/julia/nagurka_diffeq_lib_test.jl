include("nagurka_diffeq_lib.jl")
include("import_test_data.jl")
include("square_wave.jl")
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


# @testset "Comparison" begin
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
# end



function kotnik_validation()
    # fig 3 in kotnik 1998, should reach ~1.35 v in 1 microsecond

    cell_radius = 10e-6
    compare_cell = tl.Cell((0.3), (80), (0.3), (80), (3e-7), (5), (cell_radius * 2), (5e-9))

    E = 1e5
    end_time = 1e-6
    
    params = initialize_membrane_parameters(compare_cell, compare_cell, end_time, false)
    params.control_function = init_step_function(params, E, 1e-9, 1e-8)
    solution = solve_response_integrator(params)

    plot_solution(solution)
end

function talele_validation()
    tal_ref_pore_density_interpolate, tal_ref_transmembrane_voltage_interpolate, no_ep_Vm_interpolate = import_talele_test_data()

    # Note: compare specify cell radius, not diameter, as 15e-6! 
    cell_radius = 15e-6
    compare_cell = tl.Cell((1.2), (80), (0.3), (80), (3e-7), (5), (cell_radius * 2), (5e-9), py"""np.array([])""")

    E = 52000
    end_time = 1e-6
    
    params = initialize_membrane_parameters(compare_cell, compare_cell, end_time, false)
    params.control_function = init_step_function(params, E, 1e-9, 1e-8)
    solution = solve_response_integrator(params)

    plot_solution(solution)

    tau_1 = tl.first_order_schwan_time_constant(compare_cell)

    analytic_output = (1 .- exp.(-(solution.t * params.T0) / tau_1))

    formatstring = "with lines dt 1 tit "
    @gp :GP2 "set multiplot layout 2,1; set grid xtics ytics; set grid;"
    # @gp :- :GP2 1 solution.t tal_ref_transmembrane_voltage_interpolate(solution.t * params.T0)
    @gp :- :GP2 1 solution.t no_ep_Vm_interpolate(solution.t * params.T0) string(formatstring,"'Talele data'")
    @gp :- :GP2 1 solution.t col(solution, ix0_h) string(formatstring,"'Kotnik formulation'")
    @gp :- :GP2 1 solution.t analytic_output string(formatstring,"'Schwan first order'")

end


@testset "Talele thesis analytic sinusoidal validation" begin

end

# @testset "Talele transmembrane potential validation" begin


#     cell_v = cell_h = py_cell_to_julia_struct(compare_cell)

#     end_time = 1e-6
#     T0 = 1e-6
#     # this is the "apparent" nondimensionalized end time to the algorithm
#     t_f = end_time / T0

#     edge_rise_time = 1e-9
#     k = 1 / (edge_rise_time / T0)
#     x0 = 0.01

#     E = 52000

#     peak = E * (cell_h.gamma / cell_h.xi)  # ((3/2) * (compare_cell.cell_diameter / 2)) # 52 kV/m * cell diameter scale factor

#     ufun(t) = logistic_curve( t, peak, k, x0)
#     d_ufun(t) = d_logistic_curve( t, peak, k, x0)
#     d_d_ufun(t) = d_d_logistic_curve( t, peak, k, x0)

    
#     pore_solution_conductivity = 0.6 
#     params = transmembrane_params(cell_v, cell_h, nothing, nothing, nothing, t_f, nothing, nothing, tl.pore_N0, tl.pore_alpha, 
#                             tl.pore_q, tl.pore_V_ep, pore_solution_conductivity, true, T0, ufun, d_ufun, d_d_ufun)
                            
#     solution = solve_response(params)

#     N_h_course = getindex.(solution.u, Int(iN_h)+1) # note: not subtracted from N0
    
#     plot_solution(solution)

#     @gp :GP2 "set multiplot layout 2,1; set grid xtics ytics; set grid;"
#     @gp :- :GP2 1 solution.t tal_ref_transmembrane_voltage_interpolate(solution.t * T0)
#     @gp :- :GP2 1 solution.t no_ep_Vm_interpolate(solution.t * T0)
#     @gp :- :GP2 1 solution.t getindex.(solution.u, Int(ix0_h)+1)
#     @gp :- :GP2 2 solution.t tal_ref_pore_density_interpolate(solution.t * T0)
    



#     compare_cell.compute_step_response(solution.t*T0)
    
#     G_m = 1.9 # membrane conductance

#     schwan_analytic_steady_state = schwan_steady_state_with_conductivities(E, G_m, compare_cell)

#     @show schwan_analytic_steady_state

#     # @gp :GP2 "set multiplot layout 2,1; set grid xtics ytics; set grid;"
#     # @gp :- :GP2 1 solution.t getindex.(solution.u, Int(ix0_h)+1)
#     # @gp :- :GP2 1 solution.t compare_cell.step_response * 52000

#     # @gp :GP2 "set multiplot layout 5,1; set grid xtics ytics; set grid;"
#     # @gp :- :GP2 1 solution.t getindex.(solution.u, Int(ialpha)+1)
#     # @gp :- :GP2 2 solution.t getindex.(solution.u, Int(ibeta)+1)
#     # @gp :- :GP2 3 solution.t getindex.(solution.u, Int(igamma)+1)
#     # @gp :- :GP2 4 solution.t getindex.(solution.u, Int(iphi)+1)
#     # @gp :- :GP2 5 solution.t getindex.(solution.u, Int(ixi)+1)
    
#     @test_broken isapprox(getindex.(solution.u, Int(ix0_h)+1)[end], schwan_analytic)
#     @test_broken isapprox(maximum(N_h_course), 3.334e13)

# end

# @testset "DeBruin validation" begin
function debruin_validation()
    # DeBruin's formulation doesn't explicitly include the permittivity & dielectric relaxation of the cell, 
    # only the membrane capacitance.
    cell_radius = 50e-6
    E = 40000.0 # 400 V/cm to V/m
    # compare_cell = tl.Cell((0.3), (0.001), (0.3), (0.001), (3e-7), (5), (cell_radius * 2), (5e-9))
    compare_cell = tl.Cell((5.0), (0.001), (0.455), (0.001), (3e-7), (5), (cell_radius * 2), (5e-9), py"""np.array([])""")
    compare_cell.pore_solution_conductivity = 1.3

    end_time = 10e-6
    
    params = initialize_membrane_parameters(compare_cell, compare_cell, end_time, false)
    params.control_function = init_step_function(params, E, 1e-9, 1e-8)
    solution_no_ep = solve_response_integrator(params)

    params = initialize_membrane_parameters(compare_cell, compare_cell, end_time, true)
    params.control_function = init_step_function(params, E, 1e-9, 1e-8)
    solution_ep = solve_response_integrator(params)

    formatstring = "with lines dt 1 tit "
    @gp :GP2 "set multiplot layout 1,1; set grid xtics ytics; set grid;"
    @gp :- :GP2 1 solution_no_ep.t col(solution_no_ep, ix0_h) string(formatstring,"'No_ep'")
    @gp :- :GP2 1 solution_ep.t col(solution_ep, ix0_h) string(formatstring,"'ep'")

    plot_solution(solution_ep)
end




# kotnik_validation()
# talele_validation()
debruin_validation()
