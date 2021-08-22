
include("pore_lib.jl")


@testset "a" begin

    t_f = 2.0
    a = [-1.95643e-3]
    b = [1.442172e-3]

    M = 1
    m = [1.0:M;]
    t = range(epsilon, stop=t_f, length=500)

    X_t0 = 0.0
    X_tf = 1.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    virus_membrane_thickness = virus.membrane_thickness# 5e-9 # overriding temporarily! TODO: fix the root cause here!
    cell_v = cell_struct(virus.alpha, virus.beta,virus.gamma,virus.phi,virus.xi,
                virus.membrane_permittivity, virus_membrane_thickness, virus.cell_diameter)

    _X_(n) = X_(n,p,a,b,m,t_f)
    _d_X_(n) = d_X_(n,p,a,b,m,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,m,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)

    _d_V_ep(v) = d_V_ep(v, 1e15, cell_v)
    _d_d_V_ep(v,v2) = d_d_V_ep(v, v2, 1e15, cell_v)

    @gp diff(_d_V_ep.(X))/(t[end]/500)
    @gp :- _d_d_V_ep.(X, d_X)

end

