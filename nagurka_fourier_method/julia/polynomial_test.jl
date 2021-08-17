
using Test, ForwardDiff, NumericalIntegration;
using Plots;
using LinearAlgebra;

epsilon = 1e-15

function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    return (X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf)
end

# https://gcalderone.github.io/Gnuplot.jl/v1.1.0/

function P_(t, p, a, b, M, t_f)
    """
    Polynomial from polynomial_system_of_equations.py
    """
    P_t0, d_P_t0, d_d_P_t0, P_tf, d_P_tf, d_d_P_tf = p

    p_0 = P_t0
    p_1 = d_P_t0
    p_2 = d_d_P_t0/2
    p_3 = (t_f^2*(-3*d_d_P_t0 + d_d_P_tf) - 4*t_f*(3*d_P_t0 + 2*d_P_tf) - 20*P_t0 + 20*P_tf)/(2*t_f^3)
    p_4 = (t_f^2*(3*d_d_P_t0 - 2*d_d_P_tf)/2 + t_f*(8*d_P_t0 + 7*d_P_tf) + 15*P_t0 - 15*P_tf)/t_f^4
    p_5 = (t_f^2*(-d_d_P_t0 + d_d_P_tf) - 6*t_f*(d_P_t0 + d_P_tf) - 12*P_t0 + 12*P_tf)/(2*t_f^5)

    return p_0 + p_1*t + p_2*t^2 + p_3*t^3 + p_4*t^4 + p_5*t^5
end

d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> P_(n, p, a, b, M, t_f), t)
d_d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> d_P_(n, p, a, b, M, t_f), t)


@testset "Plot BCs" begin

"""
Status: works when begin and end bcs are the same.
"""
    t_f = 3.0

    M = 1
    t = range(epsilon, stop=t_f, length=60)

    P_t0 = 1.0
    P_tf = 1.0
    d_P_t0 = 4
    d_P_tf = 1.0
    d_d_P_t0 = 4.0
    d_d_P_tf = 1.0
    a = []
    b = []
    p = (P_t0, d_P_t0, d_d_P_t0, P_tf, d_P_tf, d_d_P_tf)

    _P_(n) = P_(n, p, a, b, M, t_f)
    _d_P_(n) = d_P_(n, p, a, b, M, t_f)
    _d_d_P_(n) = d_d_P_(n, p, a, b, M, t_f)

    @test isapprox(_P_(epsilon),P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_P_(t_f),P_tf, atol=1e-8, rtol=1e-8)

    @test isapprox(_d_P_(epsilon),d_P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_d_P_(t_f),d_P_tf, atol=1e-8, rtol=1e-8)

    @test isapprox(_d_d_P_(epsilon),d_d_P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_d_d_P_(t_f),d_d_P_tf, atol=1e-8, rtol=1e-8)

    plot(t, _P_.(t), show=true)

end
