
using Test;



epsilon = 1e-15


# Nagurka M, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.





function P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
    m = 1:M
    @show m
    @show sum(m.*a)

    tau = t_f
    p_0 = X_t0 - sum(a)
    p_f = X_tf - sum(a)

    return p0, pf, d_p0, d_pf, d_d_p0, d_d_pf
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

        p_0, p_f = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)
        @test isapprox(p_0, X_t0, atol=1e-8, rtol=1e-8)
        @test isapprox(p_f, X_tf, atol=1e-8, rtol=1e-8)


    end;

end;

function L_(t, a, b)

end



# function P_(t)
#
#     P = (-6 * (p_0 - p_f)
#     return
# end





function Kirk_A_optimal(t)
    return 7.289*t -6.103 + 6.696*exp(-t)-0.593*exp(t)
end

def d_Kirk_A_optimal(t):
    return 7.289 -6.696*exp(-t)-0.593*exp(t)

def d_d_Kirk_A_optimal(t):
    # not discussed in the paper - just a guess
    return 6.696*exp(-t)-0.593*exp(t)


@testset "Kirk_A_example" begin
    """
    Equation 32, Nagurka&Yen 1990, Case A
    """

    t_f = 2.0
    a = np.array([-1.95643e-3])
    b = np.array([1.442172e-3])

    M = a.shape[0]

    t = range(epsilon, t_f, 60)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    p = P_coefficients(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, M)

    X = X_(t,p,a,b,M,t_f)
    d_X = d_X_(t,p,a,b,M,t_f)
    d_d_X = d_d_X_(t,p,a,b,M,t_f)

    U = X + d_d_X

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

    assert np.allclose(X,Kirk_A_optimal(t),rtol=0.004)
    assert np.allclose(d_X,d_Kirk_A_optimal(t),rtol=0.004)
    assert np.allclose(d_d_X,d_d_Kirk_A_optimal(t),rtol=0.004)

    J = 0.5*integrate.simpson(U**2.0,t)
    assert J == pytest.approx(1.675e1)
end
