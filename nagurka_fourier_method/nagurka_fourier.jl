
using Test, ForwardDiff, NumericalIntegration;
using Plots;
using LinearAlgebra;

epsilon = 1e-15


# Nagurka M, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.


#const fdcache = ForwardDiffCache()

#maybe try sin



function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    return (X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf)
end


# function L_(t, a, b, M, t_f)
#     m = 1:M
#     return sum(a.*cos.(2*pi*(m*t)/t_f)) + sum(b.*sin.(2*pi*(m*t)/t_f))
# end;


d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, M, t_f), t)
d_d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, M, t_f), t)



d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> P_(n, p, a, b, M, t_f), t)
d_d_P_(t, p, a, b, M, t_f) = ForwardDiff.derivative(n -> d_P_(n, p, a, b, M, t_f), t)


function X_(t,p,a,b,M,t_f)
    return P_(t, p, a, b, M, t_f) + L_(t, a, b, M, t_f)
end

function d_X_(t,p,a,b,M,t_f)
    return d_P_(t, p, a, b, M, t_f) + d_L_(t, a, b, M, t_f)
end

function d_d_X_(t,p,a,b,M,t_f)
    return d_d_P_(t, p, a, b, M, t_f) + d_d_L_(t, a, b, M, t_f)
end

#
# @testset "L_" begin
#
#     t_f = 2.0
#
#     a = [1.0, -2]
#     b = [-1.0, 3]
#     M = 2
#
#     t = range(epsilon,stop=t_f,length=100)
#
#     v = (2*pi/t_f)
#
#     t_ = epsilon
#     @test isapprox(L_(t_, a, b, M, t_f), 1.0*cos(1*v*t_) + -2.0*cos(2*v*t_) + -1.0*sin(1*t_) + 3.0*sin(2*t_), atol=1e-8, rtol=1e-8)
#     t_ = 1.0
#     @test isapprox(L_(t_, a, b, M, t_f), 1.0*cos(1*t_*v) + -2.0*cos(2*t_*v) + -1.0*sin(1*t_*v) + 3.0*sin(2*t_*v), atol=1e-8, rtol=1e-8)
#
#
#     d_L_compare(t) = -1.0*v*sin(v*t) + 2.0*v*2.0*sin(2.0*t*v) + -1.0*v*cos(t*v) + 3.0*v*2.0*cos(2.0*t*v)
#     t_ = 0.5
#     @test isapprox(d_L_(t_, a, b, M, t_f), d_L_compare(t_), atol=1e-8, rtol=1e-8)
#     t_ = 2.0
#     @test isapprox(d_L_(t_, a, b, M, t_f), d_L_compare(t_), atol=1e-8, rtol=1e-8)
#
#     #
#     t_ = [0.5, 2.0]
#     d_L(t) = d_L_(t, a, b, M, t_f)
#     @test isapprox(d_L.(t_), d_L_compare.(t_), atol=1e-8, rtol=1e-8)
#
#
# end;


function Kirk_A_optimal(t)
    return 7.289.*t.-6.103.+6.696*exp.(-t).-0.593*exp.(t)
end
#
function d_Kirk_A_optimal(t)
    return 7.289.-6.696.*exp.(-t).-0.593.*exp.(t)
end

function d_d_Kirk_A_optimal(t)
    # not discussed in the paper - just a guess
    return 6.696.*exp.(-t).-0.593.*exp.(t)
end


# @testset "Kirk_A_example" begin
#     """
#     Equation 32, Nagurka&Yen 1990, Case A
#
#     no, making m=2 doesn't work
#     """
#
#     t_f = 2.0
#     a = [-1.95643e-3]
#     b = [1.442172e-3]
#
#     M = 0
#     t = range(epsilon, stop=t_f, length=60)
#
#     X_t0 = 0.0
#     X_tf = 5.0
#     d_X_t0 = 0.0
#     d_X_tf = 2.0
#     d_d_X_t0 = 6.1025137
#     d_d_X_tf = -3.4798053
#
#
#     p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#
#     _X_(n) = X_(n,p,a,b,M,t_f)
#     _d_X_(n) = d_X_(n,p,a,b,M,t_f)
#     _d_d_X_(n) = d_d_X_(n,p,a,b,M,t_f)
#     X = _X_.(t)
#     d_X = _d_X_.(t)
#     d_d_X = _d_d_X_.(t)
#
#     U = X + d_d_X
#
#     _P_(n) = P_(n, p, a, b, M, t_f)
#     plot(t,X, show=true)
#     plot!(t,_P_.(t), show=true)
#
#     plot!(t,Kirk_A_optimal(t), show=true)
#     plot!(t,d_X, show=true)
#     plot!(t,d_Kirk_A_optimal(t), show=true)
#
#     plot!(t,d_d_X, show=true)
#     plot!(t,d_d_Kirk_A_optimal(t), show=true)
#
#     @test isapprox(X,Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#     @test isapprox(d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#     @test isapprox(d_d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#
#     plot!(t,U, show=true)
#
#     J = 0.5*integrate(t, U.*U, SimpsonEven())
#     @test isapprox(J, 1.675e1, atol=1e-8, rtol=1e-8)
# end

# @testset "Kirk_C_example" begin
#     """
#     Equation 32, Nagurka&Yen 1990, Case C
#     """
#
#     t_f = 2.0
#     a = [4.3833045e-4]
#     b = [1.4015869e-3]
#
#     M = 1
#     t = range(epsilon, stop=t_f, length=60)
#
#     X_t0 = 0.0
#     X_tf = 2.3530569
#     d_X_t0 = 0.0
#     d_X_tf = 2.5293861
#     d_d_X_t0 = 1.3739699
#     d_d_X_tf = 1.9403533
#
#
#     p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#
#     _X_(n) = X_(n,p,a,b,M,t_f)
#     _d_X_(n) = d_X_(n,p,a,b,M,t_f)
#     _d_d_X_(n) = d_d_X_(n,p,a,b,M,t_f)
#     X = _X_.(t)
#     d_X = _d_X_.(t)
#     d_d_X = _d_d_X_.(t)
#
#     U = X + d_d_X
#
#     _P_(n) = P_(n, p, a, b, M, t_f)
#     _L_(n) = L_(n, a, b, M, t_f)
#
#     plot(t,X, show=true)
#     plot!(t,_P_.(t), show=true)
#     plot!(t,_L_.(t), show=true)
#
#     # plot!(t,Kirk_A_optimal(t), show=true)
#     plot!(t,d_X, show=true)
#     # plot!(t,d_Kirk_A_optimal(t), show=true)
#
#     plot!(t,d_d_X, show=true)
#     # plot!(t,d_d_Kirk_A_optimal(t), show=true)
#     #
#     # @test isapprox(X,Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#     # @test isapprox(d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#     # @test isapprox(d_d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
#
#     plot!(t,U, show=true)
#
#     J = 0.5*integrate(t, U.*U, SimpsonEven())
#     @test isapprox(J, 6.708092, atol=1e-8, rtol=1e-8)
# end


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
