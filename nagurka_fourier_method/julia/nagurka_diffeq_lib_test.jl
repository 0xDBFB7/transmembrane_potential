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



@testset "Talele validation" begin

    talale_cell = tl.Cell((1.2), (80), (0.3), (80), (3e-7), (5), (15e-6), (5e-9), py"""np.array([])""")


    end_time = 5e-6

    T0 = 1e-6 

    # this is the "apparent" nondimensionalized end time to the algorithm
    t_f = end_time / T0
    


    # edge_rise_time = 1e-7



    k = t_f * 10
    x0 = t_f / 4
    peak = 1.5

    ufun(t) = logistic_curve( t, peak, k, x0)
    d_ufun(t) = d_logistic_curve( t, peak, k, x0)
    d_d_ufun(t) = d_d_logistic_curve( t, peak, k, x0)

    params = transmembrane_params(cell_v, cell_h, a, b, p, t_f, M, m, tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep, T0, ufun, d_ufun, d_d_ufun)
    tspan = (epsilon, t_f)
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)
    solution = solve(prob, Tsit5(), dtmax = t_f / 300, maxiters= 100000, dtmin=1e-20, progress = true, progress_steps = 100)
    
    
end

