
using Test, ForwardDiff, NumericalIntegration;
using Plots;


epsilon = 1e-15


# Nagurka M, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.


#const fdcache = ForwardDiffCache()

#maybe try sin


function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    m = 1:M

    tau = t_f
    p0 = X_t0 - sum(a)
    pf = X_tf - sum(a)

    #originally + here!
    d_p0 = d_X_t0 - 2.0*sum(m.*((b.*pi)/t_f)) #reversed  sign fixed polynomial BC test?!?
    d_pf = d_X_tf - 2.0*sum(m.*((b.*pi)/t_f)) # what on earth

    d_d_p0 = d_d_X_t0 + 4.0 * sum(a.*(m.^2).*((pi^2) / (t_f^2)))
    d_d_pf = d_d_X_tf + 4.0 * sum(a.*(m.^2).*((pi^2) / (t_f^2)))

    return (p0, pf, d_p0, d_pf, d_d_p0, d_d_pf)
end

function L_(t, a, b, M, t_f)
    m = 1:M
    return sum(a.*cos.(2*pi*(m*t)/t_f)) + sum(b.*sin.(2*pi*(m*t)/t_f))
end;

d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, M, t_f), t)
d_d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, M, t_f), t)

function P_(t, p, a, b, M, t_f)
    p0, pf, d_p0, d_pf, d_d_p0, d_d_pf = p

    tau = t/t_f

    P = (-6.0*(p0 - pf) - 3.0*(d_p0 + d_pf)*t_f - (0.5*(d_d_p0 - d_d_pf)*(t_f^2)))*(tau^5)
    P += (15*(p0 - pf) + (8*d_p0 + 7*d_pf)*t_f + ((3/2.0)*d_d_p0 - d_d_pf)*(t_f^2))*(tau^4)
    P += (-10*(p0 - pf) - (6.0*d_p0 + 4.0*d_pf)*t_f - 0.5*(3.0*d_d_p0 - d_d_pf)*(t_f^2))*(tau^3)
    P += (0.5*d_d_p0*(t_f^2))*(tau^2) + (d_p0*t_f)*tau + p0

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


@testset "P_coefficients" begin

    t_f = 2.0

    a = [1.0]
    b = [-2.0]
    M = 1

    t = range(epsilon,stop=t_f,length=100)

    X_t0 = 1.0 # is it possible that this is overconstrained?
    X_tf = 2.0
    d_X_t0 = -2.0
    d_X_tf = 3.0
    d_d_X_t0 = 6
    d_d_X_tf = -3.0

    # @testset "Coefficient BCs" begin
    #     p0, pf = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    #     @test isapprox(p0, X_t0, atol=1e-8, rtol=1e-8)
    #     @test isapprox(pf, X_tf, atol=1e-8, rtol=1e-8)
    # end;

    @testset "Polynomial BCs" begin

        p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
        P(n) = P_(n, p, a, b, M, t_f)
        d_P(n) = d_P_(n, p, a, b, M, t_f)
        d_d_P(n) = d_d_P_(n, p, a, b, M, t_f)


        @test isapprox(P(epsilon), X_t0 - L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
        @test isapprox(P(t_f), X_tf - L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)

        # this appears to be what's broken.
        @test isapprox(d_P(epsilon), d_X_t0 - d_L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
        @test isapprox(d_P(t_f), d_X_tf - d_L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)

        @test isapprox(d_d_P(epsilon), d_d_X_t0 - d_d_L_(epsilon,a,b,M,t_f), atol=1e-8, rtol=1e-8)
        @test isapprox(d_d_P(t_f), d_d_X_tf - d_d_L_(t_f,a,b,M,t_f), atol=1e-8, rtol=1e-8)


    end

end;

@testset "L_" begin

    t_f = 2.0

    a = [1.0, -2]
    b = [-1.0, 3]
    M = 2

    t = range(epsilon,stop=t_f,length=100)

    v = (2*pi/t_f)

    t_ = epsilon
    @test isapprox(L_(t_, a, b, M, t_f), 1.0*cos(1*v*t_) + -2.0*cos(2*v*t_) + -1.0*sin(1*t_) + 3.0*sin(2*t_), atol=1e-8, rtol=1e-8)
    t_ = 1.0
    @test isapprox(L_(t_, a, b, M, t_f), 1.0*cos(1*t_*v) + -2.0*cos(2*t_*v) + -1.0*sin(1*t_*v) + 3.0*sin(2*t_*v), atol=1e-8, rtol=1e-8)


    d_L_compare(t) = -1.0*v*sin(v*t) + 2.0*v*2.0*sin(2.0*t*v) + -1.0*v*cos(t*v) + 3.0*v*2.0*cos(2.0*t*v)
    t_ = 0.5
    @test isapprox(d_L_(t_, a, b, M, t_f), d_L_compare(t_), atol=1e-8, rtol=1e-8)
    t_ = 2.0
    @test isapprox(d_L_(t_, a, b, M, t_f), d_L_compare(t_), atol=1e-8, rtol=1e-8)

    #
    t_ = [0.5, 2.0]
    d_L(t) = d_L_(t, a, b, M, t_f)
    @test isapprox(d_L.(t_), d_L_compare.(t_), atol=1e-8, rtol=1e-8)


end;


function Kirk_A_optimal(t)
    return 7.289.*t.-6.103.+6.696*exp.(-t).-0.593*exp.(t)
end
#
function d_Kirk_A_optimal(t)
    return 7.289 -6.696*exp.(-t)-0.593*exp.(t)
end

function d_d_Kirk_A_optimal(t)
    # not discussed in the paper - just a guess
    return 6.696*exp.(-t)-0.593*exp.(t)
end


@testset "Kirk_A_example" begin
    """
    Equation 32, Nagurka&Yen 1990, Case A
    """

    t_f = 2.0
    a = [-1.95643e-3]
    b = [1.442172e-3]

    M = 1
    t = range(epsilon, stop=t_f, length=60)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    _X_(n) = X_(n,p,a,b,M,t_f)
    _d_X_(n) = d_X_(n,p,a,b,M,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,M,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)

    U = X + d_d_X

    plot(t,X, show=true)

    plot(t,U, show=true)

    # all of these are perfect except U, need asserts
    # plt.plot(t,X)
    # plt.plot(t,P_(t, p, M))
    # plt.plot(t,
    # plt.plot(t[:-1],np.diff(X)/(1/30))
    # plt.plot(t,d_X)
    # plt.show()
    # plt.plot(t[:-2],np.diff(np.diff(X))/((1/30)**2.0))
    # plt.plot(t,d_d_X)
    # plt.plot(t,Kirk_A_optimal(t))
    # plt.plot(t,d_Kirk_A_optimal(t))
    # plt.plot(t,d_d_Kirk_A_optimal(t))
    # plt.plot(t,U)
    # plt.show()

    # this test fails. the derivatives all work out perfectly, it's just the
    # U control has a completely different shape as fig 1a 1990.
    # must be misunderstanding the U completely somehow.

    @test isapprox(X,Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
    @test isapprox(d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
    @test isapprox(d_d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)



    J = 0.5*integrate(t, U.*U, SimpsonEven())
    @test isapprox(J, 1.675e1, atol=1e-8, rtol=1e-8)
end
