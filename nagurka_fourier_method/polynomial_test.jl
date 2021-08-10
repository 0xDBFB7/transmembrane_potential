
using Test, ForwardDiff, NumericalIntegration;
using Plots;
using LinearAlgebra;

epsilon = 1e-15

function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    return (X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf)
end

function P_(t, p, a, b, M, t_f)



end

d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> P_(n, p, a, b, M, t_f), t)
d_d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> d_P_(n, p, a, b, M, t_f), t)


@testset "Plot BCs" begin

"""
Status: works when begin and end bcs are the same.
"""
    t_f = 1.0
    a = [0.0]
    b = [0.0]

    M = 1
    t = range(epsilon, stop=t_f, length=60)

    X_t0 = 1.0
    X_tf = 1.0
    d_X_t0 = 4
    d_X_tf = 1.0
    d_d_X_t0 = 4.0
    d_d_X_tf = 1.0


    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    _X_(n) = X_(n,p,a,b,M,t_f)
    _d_X_(n) = d_X_(n,p,a,b,M,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,M,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)

    @show X[end]
    @show d_X[end]
    @show d_d_X[end]

    _P_(n) = P_(n, p, a, b, M, t_f)
    _L_(n) = L_(n, a, b, M, t_f)

    plot!(t,_P_.(t), show=true)
    plot!(t,_L_.(t), show=true)
    plot(t,X, show=true)
    plot!(t,d_X, show=true)
    plot!(t,d_d_X, show=true)

end
