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


Vsc shortcuts:
tree: ctl 0
editor: ctl 1
term: ctl k
close tab: ctl w
=#

@testset "Comparison" begin
    disable_timer!(to)
    M = 3
    m = [1.0:M;]

    t_f = 1e-6

    initial_state_variables = zeros(length(instances(svars)))
    initial_state_variables[iN_v] = tl.pore_N0
    initial_state_variables[iN_h] = tl.pore_N0

    a = [0.0, 0.0, 1.0e7]
    b = [0.0, 0.0, 0.1e7]
    # a = rand(M)
    # b = rand(M)


    c = zeros(6)



    
    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = c

    X_t0 = 0.0 # override
    X_tf = 0.0
    d_X_t0 = 0.0

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    O = zeros((2*M)+6)

    virus_membrane_thickness = 5e-9 # overriding temporarily!

    cell_v = cell_struct(virus.alpha, virus.beta,virus.gamma,virus.phi,virus.xi,
                virus.membrane_permittivity, virus_membrane_thickness, virus.cell_diameter)
    cell_h = cell_struct(host_cell.alpha, host_cell.beta,host_cell.gamma,host_cell.phi,host_cell.xi,
                                    host_cell.membrane_permittivity, host_cell.membrane_thickness, host_cell.cell_diameter)
    params = transmembrane_params(cell_v, cell_h, a, b, p, t_f, M, m, tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep)
    # params.O[1:M] = a
    # params.O[M+1:(2*M)] = b
    # params.O[(2*M)+1:(2*M)+6] = c

    # params.a = O[1:M]
    # params.b = O[M+1:(2*M)]
    # params.c = O[(2*M)+1:(2*M)+6]



    tspan = (epsilon, 1e-6)
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)

    solution = solve(prob)#, RadauIIA5())
    
    solve_time = @elapsed solve(prob)#,  RadauIIA5())
    @show solve_time
    @show (solve_time / length(solution.t)) * 1e6

    integrate(solution.t, getindex.(solution.u, Int(iu0)+1))
    @btime integrate($solution.t, getindex.($solution.u, Int(iu0)+1))

    N_v_course = getindex.(solution.u, Int(iN_v)+1) .- tl.pore_N0
    N_h_course = getindex.(solution.u, Int(iN_h)+1) .- tl.pore_N0
    # Gnuplot.quitall() - plays poorly with multiplot
    @gp "set multiplot layout 3,2"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) "with lines tit 'u0'"
    @gp :- 2 solution.t getindex.(solution.u, Int(ix0_v)+1) "with lines tit 'x0_v'"
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_h)+1) "with lines tit 'x0_h'"
    @gp :- 4 solution.t N_v_course "with lines tit 'N_v'"
    @gp :- 5 solution.t N_h_course "with lines tit 'N_h'"
    # @gp :- 2 solution.t getindex.(solution.u, Int(ix1_v)+1) "with lines tit 'x0_2'"
    # @gp :- 2 solution.t getindex.(solution.u, Int(iN_h)+1) "with lines tit 'N_v'"

    
    # @gp solution.t getindex.(solution.u, Int(iu1)+1) "with lines tit 'u1'"
    # @gp solution.t getindex.(solution.u, Int(iu2)+1) "with lines tit 'u2'"
    
    # @gp solution.t getindex.(solution.u, Int(iu0)+1) "with lines tit 'N_h'"

    @testset "time" begin
        d = zeros(length(instances(svars)))
        s = zeros(length(instances(svars)))
        t=epsilon
        
        @btime transmembrane_diffeq($d,$s,$params,t2) setup=(t2 = rand())
    end

    
    

end