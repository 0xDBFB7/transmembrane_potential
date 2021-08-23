include("nagurka_diffeq_lib.jl")

using Optim
using ForwardDiff

# uncertainty stuff discussed in the unreasonable effectiveness video might be useful

const M = 20



function evaluate_control(O) 
    basetype(n) = Double64(n)

    m = [1.0:M;]
    t_f = 1e-6

    a = O[1:M]*1e4
    b = O[M+1:(2*M)]*1e4

    #Double64.
    initial_state_variables = Double64.(zeros(length(instances(svars)))) #convert(Array{BigFloat},zeros(length(instances(svars))))
    initial_state_variables[iN_v] = tl.pore_N0
    initial_state_variables[iN_h] = tl.pore_N0


    c = O[(2*M)+1:(2*M)+6]

    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = c

    X_t0 = 0.0 # override
    X_tf = 0.0
    d_X_t0 = 0.0
    d_X_t0 /= t_f
    # d_X_tf /= t_f
    d_d_X_t0 /= (t_f)
    d_d_X_tf /= (t_f)

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)
    
    virus_membrane_thickness = virus.membrane_thickness# 5e-9 # overriding temporarily! FIXME TODO

    cell_v = cell_struct(virus.alpha, virus.beta,virus.gamma,virus.phi,virus.xi,
                virus.membrane_permittivity, virus_membrane_thickness, virus.cell_diameter)
    cell_h = cell_struct(host_cell.alpha, host_cell.beta,host_cell.gamma,host_cell.phi,host_cell.xi,
                                    host_cell.membrane_permittivity, host_cell.membrane_thickness, host_cell.cell_diameter)
                        
    params = transmembrane_params(cell_v, cell_h, a, b, p, t_f, M, m, tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep)
    
    tspan = (Double64(epsilon), Double64(t_f))
    # tspan = convert.(eltype(O),tspan)
    # tspan = (epsilon, (t_f))
    prob = ODEProblem(transmembrane_diffeq,initial_state_variables,tspan,params)
    prob = remake(prob;u0=convert.(eltype(O),prob.u0))

    solution = solve(prob, Tsit5(), atol=1e-12, dtmax = t_f / 1000, progress = true, progress_steps = 500)

    N_v_course = getindex.(solution.u, Int(iN_v)+1) .- tl.pore_N0
    N_h_course = getindex.(solution.u, Int(iN_h)+1) .- tl.pore_N0

    N_v_integral = integrate(solution.t, N_v_course)
    N_h_integral = integrate(solution.t, N_h_course)

    return solution, N_v_integral, N_h_integral
end 

function cost_function(O) 
    _, N_v_integral, N_h_integral = evaluate_control(O)
    cost = -log(N_v_integral) + log(N_h_integral)
    @show ForwardDiff.value(N_v_integral), ForwardDiff.value(N_h_integral), ForwardDiff.value(cost)
    return cost
end

function optimize_coefficients() 
    O0 = Double64.(ones((2*M)+6))
    res = optimize(cost_function, O0, NelderMead())#BFGS(); autodiff = :forward)
    Ostar = Optim.minimizer(res)

    solution, _, _ = evaluate_control(Ostar)

    formatstring = "with lines ls -1 dt 1 tit "
    @gp "set multiplot layout 3,2; set grid xtics ytics; set grid;"
    @gp :- 1 solution.t getindex.(solution.u, Int(iu0)+1) string(formatstring,"'u0'")
    @gp :- 2 solution.t getindex.(solution.u, Int(iu1)+1) string(formatstring,"'u1'")
    @gp :- 3 solution.t getindex.(solution.u, Int(ix0_v)+1) string(formatstring,"'x0_v'")
    @gp :- 4 solution.t getindex.(solution.u, Int(ix0_h)+1) string(formatstring,"'x0_h'")
    @gp :- 5 solution.t N_v_course string(formatstring,"'N_v'")
    @gp :- 6 solution.t N_h_course string(formatstring,"'N_h'")
end

optimize_coefficients()

