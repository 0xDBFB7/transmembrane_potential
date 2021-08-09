
using Test, ForwardDiff;



epsilon = 1e-15


# Nagurka M, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.


#const fdcache = ForwardDiffCache()



function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    m = 1:M
    @show m
    @show sum(m.*a)

    tau = t_f
    p0 = X_t0 - sum(a)
    pf = X_tf - sum(a)

    return p0, pf #d_p0, d_pf, d_d_p0, d_d_pf
end

@testset "Setup" begin

    t_f = 2.0

    a = [-1.95643e-3]
    b = [1.442172e-3]
    M = 1

    t = range(epsilon,stop=t_f,length=100)

    X_t0 = 1.0 # is it possible that this is overconstrained?
    X_tf = 2.0
    d_X_t0 = 0
    d_X_tf = 3.0
    d_d_X_t0 = 0
    d_d_X_tf = -3.0


    @testset "Coefficient BCs" begin

        p0, pf = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
        @test isapprox(p0, X_t0, atol=1e-8, rtol=1e-8)
        @test isapprox(pf, X_tf, atol=1e-8, rtol=1e-8)


    end;

end;

function L_(t, a, b, M, t_f)
    m = 1:M
    return sum(a.*cos.(2*pi*(m*t)/t_f)) + sum(b.*sin.(2*pi*(m*t)/t_f))
end;

@testset "L_" begin

    t_f = 2.0

    a = [1.0, -2]
    b = [-1.0, 3]
    M = 2

    t = range(epsilon,stop=t_f,length=100)

    d_L_(t, a, b, M, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, M, t_f), t)

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
    correct = -1.0*v*sin.(v*t_) + 2.0*v*2.0*sin.(2.0*t_*v) + -1.0*v*cos.(t_*v) + 3.0*v*2.0*cos.(2.0*t_*v)
    @show correct
    d_L(t) = d_L_(t, a, b, M, t_f)
    @test isapprox(d_L.(t_), correct, atol=1e-8, rtol=1e-8)


end;




# function P_(t)
#
#     P = (-6 * (p_0 - p_f)
#     return
# end


#
#
#
# function Kirk_A_optimal(t)
#     return 7.289*t -6.103 + 6.696*exp(-t)-0.593*exp(t)
# end
#
# def d_Kirk_A_optimal(t):
#     return 7.289 -6.696*exp(-t)-0.593*exp(t)
#
# def d_d_Kirk_A_optimal(t):
#     # not discussed in the paper - just a guess
#     return 6.696*exp(-t)-0.593*exp(t)
#
#
# @testset "Kirk_A_example" begin
#     """
#     Equation 32, Nagurka&Yen 1990, Case A
#     """
#
#     t_f = 2.0
#     a = np.array([-1.95643e-3])
#     b = np.array([1.442172e-3])
#
#     M = a.shape[0]
#
#     t = range(epsilon, t_f, 60)
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
#     X = X_(t,p,a,b,M,t_f)
#     d_X = d_X_(t,p,a,b,M,t_f)
#     d_d_X = d_d_X_(t,p,a,b,M,t_f)
#
#     U = X + d_d_X
#
#     # all of these are perfect except U, need asserts
#     # plt.plot(t,X)
#     # plt.plot(t,P_(t, p, M))
#     # plt.plot(t,
#     # plt.plot(t[:-1],np.diff(X)/(1/30))
#     # plt.plot(t,d_X)
#     # plt.show()
#     # plt.plot(t[:-2],np.diff(np.diff(X))/((1/30)**2.0))
#     # plt.plot(t,d_d_X)
#     # plt.plot(t,Kirk_A_optimal(t))
#     # plt.plot(t,d_Kirk_A_optimal(t))
#     # plt.plot(t,d_d_Kirk_A_optimal(t))
#     # plt.plot(t,U)
#     # plt.show()
#
#     # this test fails. the derivatives all work out perfectly, it's just the
#     # U control has a completely different shape as fig 1a 1990.
#     # must be misunderstanding the U completely somehow.
#
#     assert np.allclose(X,Kirk_A_optimal(t),rtol=0.004)
#     assert np.allclose(d_X,d_Kirk_A_optimal(t),rtol=0.004)
#     assert np.allclose(d_d_X,d_d_Kirk_A_optimal(t),rtol=0.004)
#
#     J = 0.5*integrate.simpson(U**2.0,t)
#     assert J == pytest.approx(1.675e1)
# end
