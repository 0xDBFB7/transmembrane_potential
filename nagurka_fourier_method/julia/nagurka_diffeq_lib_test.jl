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

# @testset "current" begin
#     V_m = 5.0
#     N = 1e20
#     cell_h = cell_struct(host_cell.alpha, host_cell.beta,host_cell.gamma,host_cell.phi,host_cell.xi,
#                                     host_cell.membrane_permittivity, host_cell.membrane_thickness, host_cell.cell_diameter)
#     @show dVdt_tm_potential_pore_current(V_m, N, cell_h)

#     vs = [0.0:3.0]
#     dVdt(k, N, cell_h) = dVdt_tm_potential_pore_current(k, N, cell_h)

#     @gp V_m / 
    
# end

# staticarrays.jl




@testset "Comparison" begin
    # M = 10
    # a = (zeros(M))
    # b = (zeros(M))
    # a[9] = 0.5e7
    # b[4] = -0.1e7
    
    M = 3
    a = [0.0, 0.0, 1.0e4]
    b = [0.0, 0.0, 0.1e4]

    # M = 10
    # a = rand(M)
    # b = rand(M)
    
    # M=30
    # a = (zeros(M))
    # b = (zeros(M))
    # for i in [1.0:2.0:M;]
    #     b[Int(i)] = (1/(i))*80000*(50/20)
    # end
    
    
    m = [1.0:M;]
    t_f = 1e-6

    basetype(n) = Double64(n)

    initial_state_variables = basetype.(zeros(length(instances(svars)))) #convert(Array{BigFloat},zeros(length(instances(svars))))
    initial_state_variables[iN_v] = tl.pore_N0
    initial_state_variables[iN_h] = tl.pore_N0


    c = zeros(6)

    
    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = c

    X_t0 = 0.0 # override
    X_tf = 0.0
    d_X_t0 = 0.0

    d_X_t0 = d_L_(epsilon, a, b, m, t_f)
    d_X_tf = d_L_(t_f, a, b, m, t_f)
    d_d_X_t0 = d_d_L_(epsilon, a, b, m, t_f)
    d_d_X_tf = d_d_L_(t_f, a, b, m, t_f)

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    O = zeros((2*M)+6)

    
    #local scope for speed
    virus_membrane_thickness = virus.membrane_thickness# 5e-9 # overriding temporarily! FIXME TODO
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

    # tspan = (basetype(epsilon), basetype(t_f/1000.0)) # precompiles everything
    # prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)

    _solve() = solve(prob, RadauIIA5(), atol=1e-8, dtmax = t_f / 1000, progress = true, progress_steps = 500)

    
    tspan = (Double64(epsilon), Double64(t_f))
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)

    # @btime solve($prob)#, RadauIIA5(), dtmax = t_f / 100) #, dtmin=1e-15, atol=1e-8, rtol=1e-8, #RadauIIA5(), dtmax = t_f / 100) # dtmin=1e-13,
    #@benchmark doesn't print off the repl?
    

    tspan = ((epsilon), (1e-6))
    solution = _solve()

    begin_t = time()

    solution = _solve()

    solve_time = (time()-begin_t)/10.0
    @show solve_time
    @show (solve_time / length(solution.t)) * 1e6
    
    N_v_course = getindex.(solution.u, Int(iN_v)+1) .- tl.pore_N0
    N_h_course = getindex.(solution.u, Int(iN_h)+1) .- tl.pore_N0
    # Gnuplot.quitall() - plays poorly with multiplot
    formatstring = "with lines ls -1 dt 1 tit "
    @gp "set multiplot layout 3,2; set grid xtics ytics; set grid;"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) string(formatstring,"'u0'")
    @gp :- 2 solution.t getindex.(solution.u, Int(iu1)+1) string(formatstring,"'u1'")
    # @gp :- 3 solution.t getindex.(solution.u, Int(ix1_v)+1) string(formatstring,"'x1_v'")
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_v)+1) string(formatstring,"'x0_v'")
    @gp :- 4 solution.t getindex.(solution.u, Int(ix0_h)+1) string(formatstring,"'x0_h'")
    @gp :- 5 solution.t N_v_course string(formatstring,"'N_v'")
    @gp :- 6 solution.t N_h_course string(formatstring,"'N_h'")
    
    #this is just confusing now because it's the second derivative
    # @gp :- 6 solution.t getindex.(solution.u, Int(iI_ep_v)+1) string(formatstring,"'I_ep_v'")
    
    #save(term="svg size 1000 1000",output="runs/both_derivatives_V_ep.svg")
    
    # @gp solution.t getindex.(solution.u, Int(iu1)+1) "with lines tit 'u1'"
    # @gp solution.t getindex.(solution.u, Int(iu2)+1) "with lines tit 'u2'"
    
    # @gp solution.t getindex.(solution.u, Int(iu0)+1) "with lines tit 'N_h'"
    
    @testset "time" begin
        d = zeros(length(instances(svars)))
        s = zeros(length(instances(svars)))
        t=epsilon
        
        single_timestep_time = @belapsed transmembrane_diffeq($d,$s,$params,t2) setup=(t2 = rand()) seconds=0.5
        @show single_timestep_time
        @show length(solution.t) * single_timestep_time
    end
    
    # integrate(solution.t, getindex.(solution.u, Int(iu0)+1))
    # @btime integrate($solution.t, getindex.($solution.u, Int(iu0)+1))
    
end
