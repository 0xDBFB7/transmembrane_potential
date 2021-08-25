include("nagurka_diffeq_lib.jl")

using Optim
using ForwardDiff
# using IntervalArithmetic, IntervalOptimisation
using BlackBoxOptim
# uncertainty stuff discussed in the unreasonable effectiveness video might be useful
using Printf

M = 30



# Nondimensionalize!
# No reason to have coefficients on the N - it's just an integral, doesn't need to explode...
# can leave out the decay term of the pore differential equation, since  b

function evaluate_control(O) 
    basetype(n) = Double64(n)

    m = [1.0:M;]
    # t_f =abs(O[end]*1e-4)
    t_f = 10e-8
    
    a = O[1:M]
    b = O[M+1:(2*M)]

    a /= a
    b /= b
    # ab = a.+b


    # a /= ab
    # b /= ab

    a *= 1e6 / M
    b *= 1e6 / M

    # divide O by a+b to maintain relative ab scaling?

    #Double64.
    initial_state_variables = Double64.(zeros(length(instances(svars)))) #convert(Array{BigFloat},zeros(length(instances(svars))))
    initial_state_variables[iN_v] = tl.pore_N0
    initial_state_variables[iN_h] = tl.pore_N0


    c = O[(2*M)+1:(2*M)+6]

    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = c

    X_t0 = 0.0 # override
    X_tf = 0.0
    d_X_t0 = 0.0
    d_X_tf = 0.0
    # d_d_X_t0 = 0.0
    # d_d_X_tf = 0.0

    # d_X_t0 /= t_f
    # d_X_tf /= t_f
    # d_d_X_t0 /= (t_f*t_f)
    # d_d_X_tf /= (t_f*t_f)

    d_X_t0 = 0.1*d_L_(epsilon, a, b, m, t_f)
    d_X_tf = 0.1*d_L_(t_f, a, b, m, t_f)
    d_d_X_t0 = d_d_L_(epsilon, a, b, m, t_f)
    d_d_X_tf = d_d_L_(t_f, a, b, m, t_f)
    
    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)
    
    virus_membrane_thickness = virus.membrane_thickness# 5e-9 # overriding temporarily! FIXME TODO

    cell_v = cell_struct(virus.alpha, virus.beta,virus.gamma,virus.phi,virus.xi,
                virus.membrane_permittivity, virus_membrane_thickness, virus.cell_diameter)
    cell_h = cell_struct(host_cell.alpha, host_cell.beta,host_cell.gamma,host_cell.phi,host_cell.xi,
                                    host_cell.membrane_permittivity, host_cell.membrane_thickness, host_cell.cell_diameter)
                        
    params = transmembrane_params(cell_v, cell_h, a, b, p, t_f, M, m, tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep)
    
    # tspan = (Double64(epsilon), Double64(t_f))
    # tspan = convert.(eltype(O),tspan)
    tspan = (epsilon, (t_f))
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)
    # prob = remake(prob;u0=convert.(eltype(O),prob.u0), tspan=tspan) # needed for Dual numbers
    # prob = remake(prob; tspan=tspan)

    solution = solve(prob, Tsit5(), atol=1e-7, dtmin=1e-20, dtmax = t_f / 300, progress = false, progress_steps = 500)
    # rk4 usually errorsout! Tsit5 seems to usually work
    # Vern9 is a bit faster on double64s.

    N_v_course = getindex.(solution.u, Int(iN_v)+1) .- tl.pore_N0
    N_h_course = getindex.(solution.u, Int(iN_h)+1) .- tl.pore_N0

    N_v_integral = integrate(solution.t, N_v_course)/t_f
    N_h_integral = integrate(solution.t, N_h_course)/t_f

    x0_v_integral = integrate(solution.t, getindex.(solution.u, Int(ix0_v)+1))/t_f
    x0_h_integral = integrate(solution.t, getindex.(solution.u, Int(ix0_h)+1))/t_f

    return solution, N_v_integral, N_h_integral, x0_v_integral, x0_h_integral
end 

function cost_function(O) 
    _, N_v_integral, N_h_integral, x0_v_integral, x0_h_integral = evaluate_control(O)
    # cost = -log(N_v_integral) + log(N_h_integral) #+ 1/N_v_integral
    cost = -log(N_v_integral) + log(N_h_integral) #+ 1/N_v_integral
    # cost = -N_v_integral + N_h_integral
    @printf("N_v: %.3e | N_h: %.3e | (%.3e) | x0_v: %.3e | x0_h: %.3e | (%.3e) | Cost: %.3e\n",ForwardDiff.value(N_v_integral), ForwardDiff.value(N_h_integral),
                                        (ForwardDiff.value(N_h_integral)/ForwardDiff.value(N_v_integral)),
                                        ForwardDiff.value(x0_v_integral), ForwardDiff.value(x0_h_integral), 
                                        (ForwardDiff.value(x0_h_integral)/ForwardDiff.value(x0_v_integral)),
                                         ForwardDiff.value(cost))
    return cost
end

function optimize_coefficients()#global_Ostar) 
    # O0 = global_Ostar
    O0 = (ones((2*M)+7)) # Double64. - needed for autodiff-based 
    O0[(2*M)+1:(2*M)+6] .= 0.0
    res = optimize(cost_function, O0,  NelderMead(), Optim.Options(iterations = 1, show_trace=false))
                                #  autodiff = :forward)
    # BFGS(); autodiff = :forward)
    Ostar = Optim.minimizer(res)

    # interval_ = IntervalBox(-1e2..1e2, (2*M)+7)
    # _, Ostar = minimise(cost_function, interval_, tol=1e-2)

    # res = bboptimize(cost_function; SearchRange = (-1.0, 1.0), NumDimensions = ((2*M)+7), method=:generating_set_search, MaxTime = 100.0)
    # need to wrap cost in float64()
    # Ostar = best_candidate(res)

    #=
    The problem with gradien tmethods without a log cost function
    is that, if one is 1e20 and the other 1e3, the gradient will be *huge*. 
    Non-autodiff CongugateGradient works really well but gets stuck at 35
    =#


    solution, _, _,_,_ = evaluate_control(Ostar)

    formatstring = "with lines ls -1 dt 1 tit "
    @gp "set multiplot layout 3,2; set grid xtics ytics; set grid;"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) string(formatstring,"'u0'")
    @gp :- 2 solution.t getindex.(solution.u, Int(iu1)+1) string(formatstring,"'u1'")
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_v)+1) string(formatstring,"'x0_v'")
    @gp :- 4 solution.t getindex.(solution.u, Int(ix0_h)+1) string(formatstring,"'x0_h'")
    @gp :- 5 solution.t (getindex.(solution.u, Int(iN_v)+1).- tl.pore_N0) string(formatstring,"'N_v'")
    @gp :- 6 solution.t (getindex.(solution.u, Int(iN_h)+1).- tl.pore_N0) string(formatstring,"'N_h'")

    return Ostar
end

optimize_coefficients()



