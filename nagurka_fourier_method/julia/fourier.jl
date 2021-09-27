include("nagurka_diffeq_lib.jl")

using Optim
using ForwardDiff
using NumericalIntegration
# using IntervalArithmetic, IntervalOptimisation
using BlackBoxOptim
# uncertainty stuff discussed in the unreasonable effectiveness video might be useful
using Printf

M = 10

# pkg activate dev_nagurka

# using dev version of OrdinaryDiffEq


# Nondimensionalize!
# No reason to have coefficients on the N - it's just an integral, doesn't need to explode...
# can leave out the decay term of the pore differential equation, since  b

function evaluate_control(O) 
    # basetype(n) = Double64(n)

    end_time = 1e-6

    a = O[1:M]
    b = O[M+1:(2*M)]
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

    
    virus.pore_solution_conductivity = 0.15
    host_cell.pore_solution_conductivity = 0.6
    
    params = initialize_membrane_parameters(virus, host_cell, end_time, true)

    m = [1.0:M;]
    # d_X_t0 = d_L_(epsilon, a, b, m, t_f)
    # d_X_tf = d_L_(t_f, a, b, m, t_f)
    d_d_X_t0 = d_d_L_(epsilon, a, b, m, params.t_f)
    d_d_X_tf = d_d_L_(params.t_f, a, b, m, params.t_f)


    params.control_function = init_fourier_parametrization(params, a, b, M, X_t0,
                                                 d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf)
    solution = solve_response_integrator(params)

    N_v_course = getindex.(solution.u, Int(iN_v)+1) .- tl.pore_N0
    N_h_course = getindex.(solution.u, Int(iN_h)+1) .- tl.pore_N0

    N_v_integral = NumericalIntegration.integrate(solution.t, N_v_course)/params.t_f
    N_h_integral = NumericalIntegration.integrate(solution.t, N_h_course)/params.t_f

    x0_v_integral = NumericalIntegration.integrate(solution.t, getindex.(solution.u, Int(ix0_v)+1))/params.t_f
    x0_h_integral = NumericalIntegration.integrate(solution.t, getindex.(solution.u, Int(ix0_h)+1))/params.t_f

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
    O0 = (ones((2*M)+7))  # Double64. - needed for autodiff-based 
    # O0[(2*M)+1:(2*M)+6] .= 0.0
    O0[M+1:(2*M)] .= 0
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
    The problem with gradient methods without a log cost function
    is that, if one is 1e20 and the other 1e3, the gradient will be *huge*. 
    Non-autodiff CongugateGradient works really well but gets stuck at 35
    =#


    solution, _, _,_,_ = evaluate_control(Ostar)

    return Ostar
end

# optimize_coefficients()

# a,b = square_wave(M)
Ostar = (ones((2*M)+7)) .* -1.0 / (2*M) * 2.0

# twiddling 

Ostar[(2*M)+1:(2*M)+6] .= 0.0
# Ostar[M+1:(2*M)] .= 0
solution, _, _,_,_ = evaluate_control(Ostar)
plot_solution(solution)
solution_ep = solution
@gp :GP3 "set multiplot layout 5,1; set grid xtics ytics; set grid;"
@gp :- :GP3 1 solution_ep.t getindex.(solution_ep.u, Int(ialpha)+1)
@gp :- :GP3 2 solution_ep.t getindex.(solution_ep.u, Int(ibeta)+1)
@gp :- :GP3 3 solution_ep.t getindex.(solution_ep.u, Int(igamma)+1)
@gp :- :GP3 4 solution_ep.t getindex.(solution_ep.u, Int(iphi)+1)
@gp :- :GP3 5 solution_ep.t getindex.(solution_ep.u, Int(ixi)+1)
@show virus.gamma_ep