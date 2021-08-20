include("nagurka_diffeq_lib.jl")

BenchmarkTools.DEFAULT_PARAMETERS.samples = 30

#=

To do mem alloc analysis, 

julia --track-allocation=user
<run function> - precompiles
Profile.clear_malloc_data()
<run function>


https://github.com/JuliaApproximation/FastGaussQuadrature.jl apparently
is better if you can chose where to sample the thing, like an odesolver lets you.

=#

@testset "Comparison" begin
    M = 3
    m = [1.0:M;]

    t_f = 1e-6

    initial_state_variables = zeros(length(instances(svars)))
    initial_state_variables[iN_v] = pore_N0
    initial_state_variables[iN_h] = pore_N0

    a = [0.0, 0.0, 1.0]
    b = [0.0, 0.0, 0.1]
    

    c = zeros(6)

    O = zeros((2*M)+6)
    params = transmembrane_params(O, t_f, M, m)
    params.O[1:M] = a
    params.O[M+1:(2*M)] = b
    params.O[(2*M)+1:(2*M)+6] = c


    tspan = (epsilon, 1e-6)
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)

    solution = solve(prob)

    @time solution = solve(prob)

    integrate(solution.t, getindex.(solution.u, Int(iu0)+1))
    @btime integrate($solution.t, getindex.($solution.u, Int(iu0)+1))

    # Gnuplot.quitall() - plays poorly with multiplot
    @gp "set multiplot layout 3,1"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) "with lines tit 'u0'"
    @gp :- 2 solution.t getindex.(solution.u, Int(ix0_v)+1) "with lines tit 'x0_v'"
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_h)+1) "with lines tit 'x0_h'"

    # @gp :- 2 solution.t getindex.(solution.u, Int(ix1_v)+1) "with lines tit 'x0_2'"
    # @gp :- 2 solution.t getindex.(solution.u, Int(iN_h)+1) "with lines tit 'N_v'"

    
    # @gp solution.t getindex.(solution.u, Int(iu1)+1) "with lines tit 'u1'"
    # @gp solution.t getindex.(solution.u, Int(iu2)+1) "with lines tit 'u2'"
    
    # @gp solution.t getindex.(solution.u, Int(iu0)+1) "with lines tit 'N_h'"

    @testset "time" begin
        d = zeros(length(instances(svars)))
        s = zeros(length(instances(svars)))
        t=epsilon
        @btime transmembrane_diffeq($d,$s,$params,$t) 
    end
end