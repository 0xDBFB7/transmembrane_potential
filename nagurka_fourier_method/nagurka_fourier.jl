
using Test, ForwardDiff, NumericalIntegration;
using Plots;
using LinearAlgebra;

epsilon = 1e-15


# Nagurka M, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.


#const fdcache = ForwardDiffCache()

#maybe try sin


# function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#     m = 1:M
#
#     tau = t/t_f
#     p0 = X_t0 - sum(a)
#     pf = X_tf - sum(a)
#
#     #originally + here!
#     d_p0 = d_X_t0 + 2.0*sum(m.*((b.*pi)/t_f)) #reversed  sign fixed polynomial BC test?!?
#     d_pf = d_X_tf + 2.0*sum(m.*((b.*pi)/t_f)) # what on earth
#
#     d_d_p0 = d_d_X_t0 + 4.0 * sum(a.*(m.^2).*((pi^2) / (t_f^2)))
#     d_d_pf = d_d_X_tf + 4.0 * sum(a.*(m.^2).*((pi^2) / (t_f^2)))
#
#     return (p0, pf, d_p0, d_pf, d_d_p0, d_d_pf)
# end

function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    return (X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf)
end


# function L_(t, a, b, M, t_f)
#     m = 1:M
#     return sum(a.*cos.(2*pi*(m*t)/t_f)) + sum(b.*sin.(2*pi*(m*t)/t_f))
# end;

function L_(t, a, b, M, t_f)
    """
    Yen V, Nagurka ML.
    Fourier-based optimal control approach for structural systems.
    Journal of Guidance, Control, and Dynamics 1990;13:265–76.
    https://doi.org/10.2514/3.20546.
    or
    Yen V, Nagurka ML.
    A Fourier-Based Optimal Control Approach for Structural Systems.
    1988 American Control Conference, 1988, p. 2082–7.
    https://doi.org/10.23919/ACC.1988.4790067.
    """
    m = [1.0:M;]
    #m = 1:M
    # alpha_k = -1 .+ 4(m.^2.0)*(pi^2)*((t/t_f)^2 - (2*(t/t_f)^3) + (t/t_f)^4) .+ cos.(2*pi*(m.*t)/t_f)
    # beta_k = 2.0*m.*pi*(-(t/t_f) + 10*(t/t_f)^3 - 15*(t/t_f)^4 + 6*(t/t_f)^5) .+ sin.(2*pi*(m.*t)/t_f)

    tau = t/t_f
    vk = 2 * m * pi
    alpha_k = -1 .+ (vk.^2)*(0.5*tau^2 - tau^3 + 0.5 * tau^4) .+ cos.(vk .* tau)
    beta_k = vk.* (- tau + 10 * tau^3 - 15*tau^4 + 6*tau^5) .+ sin.(vk .* tau)

    return sum(a.*alpha_k) + sum(b.*beta_k)

end;

d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, M, t_f), t)
d_d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, M, t_f), t)

# function P_(t, p, a, b, M, t_f)
#     p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = p
#
#     tau = t/t_f
#
#     P = (-6.0*(p0 - pf) - 3.0*(d_p0 + d_pf)*t_f - (0.5*(d_d_p0 - d_d_pf)*(t_f^2)))*(tau^5)
#     P += (15*(p0 - pf) + (8*d_p0 + 7*d_pf)*t_f + ((3/2.0)*d_d_p0 - d_d_pf)*(t_f^2))*(tau^4)
#     P += (-10*(p0 - pf) - (6.0*d_p0 + 4.0*d_pf)*t_f - 0.5*(3.0*d_d_p0 - d_d_pf)*(t_f^2))*(tau^3)
#     P += (0.5*d_d_p0*(t_f^2))*(tau^2) + (d_p0*t_f)*tau + p0
#
#     return P
# end

function P_(t, p, a, b, M, t_f)
    """

    Yen V, Nagurka ML.
    Fourier-based optimal control approach for structural systems.
    Journal of Guidance, Control, and Dynamics 1990;13:265–76.
    https://doi.org/10.2514/3.20546.


    THERE IS A TYPO in the equivalent p = p0 .... equation in
    Yen V, Nagurka ML.
    A Fourier-Based Optimal Control Approach for Structural Systems.
    1988 American Control Conference, 1988, p. 2082–7.
    https://doi.org/10.23919/ACC.1988.4790067.
    specifically the p1 term should be multiplied by d_d_X_t0

    Three-poly in
    Yen V, Nagurka ML.
    Multiple-segment fourier-based approach for linear quadratic optimal control 1989.
    https://doi.org/10.1109/ICCON.1989.770669.

    """
    X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf = p

    tau = t/t_f

    p0 = X_t0 + X_t0*t + (-10*X_t0 - 6*d_X_t0*t_f)*(tau)^3
            + (15*X_t0 + d_X_t0*t_f)*(tau)^4 + (-6*X_t0 - 3*d_X_t0*t_f)*(tau)^5

    p1 = 0.5*(t_f^2)*((tau^2) - 3*tau^3 + 3*tau^4 - tau^5)
    p2 = (10*tau^3 - 15*tau^4 + 6*tau^5)
    p3 = t_f*(-4*tau^3 + 7*tau^4 - 3*tau^5)
    p4 = 0.5*(t_f^2)*(tau^3 - 2*tau^4 + tau^5)

    P = p0 + p1*d_d_X_t0 + p2*X_tf + p3*d_X_tf + p4 * d_d_X_tf

    return P
end

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
# @testset "P_coefficients" begin
#
#     t_f = 2.0
#
#     a = [1.0]
#     b = [-2.0]
#     M = 1
#
#     t = range(epsilon,stop=t_f,length=100)
#
#     X_t0 = 1.0 # is it possible that this is overconstrained?
#     X_tf = 2.0
#     d_X_t0 = -2.0
#     d_X_tf = 3.0
#     d_d_X_t0 = 6
#     d_d_X_tf = -3.0
#
#     # @testset "Coefficient BCs" begin
#     #     p0, pf = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#     #     @test isapprox(p0, X_t0, atol=1e-8, rtol=1e-8)
#     #     @test isapprox(pf, X_tf, atol=1e-8, rtol=1e-8)
#     # end;
#
#     @testset "Polynomial BCs" begin
#
#         p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
#         P(n) = P_(n, p, a, b, M, t_f)
#         d_P(n) = d_P_(n, p, a, b, M, t_f)
#         d_d_P(n) = d_d_P_(n, p, a, b, M, t_f)
#
#
#         @test isapprox(P(epsilon), X_t0 - L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#         @test isapprox(P(t_f), X_tf - L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#
#         # this appears to be what's broken.
#         @test isapprox(d_P(epsilon), d_X_t0 - d_L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#         @test isapprox(d_P(t_f), d_X_tf - d_L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#
#         @test isapprox(d_d_P(epsilon), d_d_X_t0 - d_d_L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#         @test isapprox(d_d_P(t_f), d_d_X_tf - d_d_L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)
#
#
#     end
#
# end;
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
