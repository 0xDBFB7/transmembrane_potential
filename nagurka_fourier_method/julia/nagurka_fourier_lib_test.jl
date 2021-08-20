include("nagurka_fourier_lib.jl")


using ForwardDiff


autodiff_d_L_(t, a, b, m, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, m, t_f), t)
autodiff_d_d_L_(t, a, b, m, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, m, t_f), t)

autodiff_d_P_(t, P_BCs, a, b, m, t_f) = ForwardDiff.derivative(n -> P_(n, P_BCs, a, b, m, t_f), t)
autodiff_d_d_P_(t, P_BCs, a, b, m, t_f) = ForwardDiff.derivative(n -> d_P_(n, P_BCs, a, b, m, t_f), t)


@testset "Polynomial BCs" begin

"""
Status: works when begin and end bcs are the same.
"""
    t_f = 3.0

    M = 1
    m = [1.0:M;]

    t = range(epsilon, stop=t_f, length=60)


    P_t0 = 1.0
    P_tf = 1.0
    d_P_t0 = 4
    d_P_tf = 1.0
    d_d_P_t0 = 4.0
    d_d_P_tf = 1.0
    a = []
    b = []
    P_BCs = (P_t0, d_P_t0, d_d_P_t0, P_tf, d_P_tf, d_d_P_tf)

    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    _P_(n) = P_(n, p, a, b, m, t_f)
    _d_P_(n) = d_P_(n, p, a, b, m, t_f)
    _d_d_P_(n) = d_d_P_(n, p, a, b, m, t_f)

    @test isapprox(_P_(epsilon),P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_P_(t_f),P_tf, atol=1e-8, rtol=1e-8)

    @test isapprox(_d_P_(epsilon),d_P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_d_P_(t_f),d_P_tf, atol=1e-8, rtol=1e-8)

    @test isapprox(_d_d_P_(epsilon),d_d_P_t0, atol=1e-8, rtol=1e-8)
    @test isapprox(_d_d_P_(t_f),d_d_P_tf, atol=1e-8, rtol=1e-8)

    # plot(t, _P_.(t), show=true)

end

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


@testset "Kirk_A_example" begin
    """
    Equation 32, Nagurka&Yen 1990, Case A

    """

    t_f = 2.0
    a = [-1.95643e-3]
    b = [1.442172e-3]

    M = 1
    m = [1.0:M;]
    t = range(epsilon, stop=t_f, length=60)

    X_t0 = 0.0
    X_tf = 5.0
    d_X_t0 = 0.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.1025137
    d_d_X_tf = -3.4798053


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)


    _X_(n) = X_(n,p,a,b,m,t_f)

    _d_X_(n) = d_X_(n,p,a,b,m,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,m,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)

    U = X + d_d_X



    @btime X_($epsilon,$p,$a,$b,$m,$t_f)

    # _P_(n) = P_(n, P_BCs, a, b, m, t_f)
    # plot(t,X, show=true)
    # plot!(t,_P_.(t), show=true)
    #
    # plot!(t,Kirk_A_optimal(t), show=true)
    # plot!(t,d_X, show=true)
    # plot!(t,d_Kirk_A_optimal(t), show=true)
    #
    # plot!(t,d_d_X, show=true)
    # plot!(t,d_d_Kirk_A_optimal(t), show=true)
    #

    # plot!(t,X-Kirk_A_optimal(t), show=true)

    @test_broken isapprox(X,Kirk_A_optimal(t), atol=(0.04/100.0), rtol=(0.04/100.0))
    @test_broken isapprox(d_X,d_Kirk_A_optimal(t), atol=(0.04/100.0), rtol=(0.04/100.0))
    @test_broken isapprox(d_d_X,d_d_Kirk_A_optimal(t), atol=(0.04/100.0), rtol=(0.04/100.0))

    # plot!(t,U, show=true)

    J = 0.5*integrate(t, U.*U, SimpsonEven())
    @show J
    @test_broken isapprox(J, 1.675e1, atol=1e-8, rtol=1e-8)
end

macro name(x)
    quote
        $(string(x))
    end
end


@testset "Kirk_C_example" begin
    """
    Equation 32, Nagurka&Yen 1990, Case C
    """

    t_f = 2.0
    a = [4.3833045e-4]
    b = [1.4015869e-3]

    M = 1
    m = [1.0:M;]
    t = range(epsilon, stop=t_f, length=60)

    X_t0 = 0.0
    X_tf = 2.3530569
    d_X_t0 = 0.0
    d_X_tf = 2.5293861
    d_d_X_t0 = 1.3739699
    d_d_X_tf = 1.9403533


    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    _X_(n) = X_(n,p,a,b,m,t_f)
    _d_X_(n) = d_X_(n,p,a,b,m,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,m,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)

    U = X + d_d_X

    _P_(n) = P_(n, p, a, b, m, t_f)
    _L_(n) = L_(n, a, b, m, t_f)

    # plot(t,X, show=true, label=@name X)
    # plot!(t,_P_.(t), show=true, label=@name  _P_.(t))
    # plot!(t,_L_.(t), show=true, label=@name  _L_.(t))

    # plot!(t,Kirk_A_optimal(t), show=true)
    # plot!(t, d_X, show=true, label=@name  d_X)
    # plot!(t,d_Kirk_A_optimal(t), show=true)
    # plot!(t, d_d_X, show=true, label=@name  d_d_X)
    # plot!(t,d_d_Kirk_A_optimal(t), show=true)
    #
    # @test isapprox(X,Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
    # @test isapprox(d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
    # @test isapprox(d_d_X,d_d_Kirk_A_optimal(t), atol=1e-8, rtol=1e-8)
    # plot!(t, U, show=true, label=@name  U)


    J = 0.5*integrate(t, U.*U, SimpsonEven())
    @test_broken isapprox(J, 6.708092, atol=1e-8, rtol=1e-8)
end


@testset "Verify BCs" begin

    t_f = 3.0
    a = [1.0, 3]
    b = [-2.0, -1]

    M = 2
    m = [1.0:M;]
    t = range(epsilon, stop=t_f, length=60)

    X_t0 = 1.0
    X_tf = 5.0
    d_X_t0 = 4.0
    d_X_tf = 2.0
    d_d_X_t0 = 6.0
    d_d_X_tf = 3.0

    P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    p = P_BCs_to_p_coefficients(P_BCs, t_f)

    _X_(n) = X_(n,p,a,b,m,t_f)
    _d_X_(n) = d_X_(n,p,a,b,m,t_f)
    _d_d_X_(n) = d_d_X_(n,p,a,b,m,t_f)
    X = _X_.(t)
    d_X = _d_X_.(t)
    d_d_X = _d_d_X_.(t)


    @testset "Autodiff" begin
        @test isapprox(X[begin], X_t0, atol=1e-8, rtol=1e-8)
        @test isapprox(X[end], X_tf, atol=1e-8, rtol=1e-8)
        @test isapprox(d_X[begin], d_X_t0, atol=1e-8, rtol=1e-8)
        @test isapprox(d_X[end], d_X_tf, atol=1e-8, rtol=1e-8)
        @test isapprox(d_d_X[begin], d_d_X_t0, atol=1e-8, rtol=1e-8)
        @test isapprox(d_d_X[end], d_d_X_tf, atol=1e-8, rtol=1e-8)
    end

    _d_P_(n) = d_P_(n, p, a, b, m, t_f)
    _d_d_P_(n) = d_d_P_(n, p, a, b, m, t_f)

    _d_L_(n) = d_L_(n, a, b, m, t_f)
    _d_d_L_(n) = d_d_L_(n, a, b, m, t_f)


    @testset "analytic" begin
        a_d_P_(n) = autodiff_d_P_(n, p, a, b, m, t_f)
        a_d_d_P_(n) = autodiff_d_d_P_(n, p, a, b, m, t_f)

        @test isapprox(a_d_P_.(t), _d_P_.(t), atol=1e-8, rtol=1e-8)
        @test isapprox(a_d_d_P_.(t), _d_d_P_.(t), atol=1e-8, rtol=1e-8)

        a_d_L_(n) = autodiff_d_L_(n, a, b, m, t_f)
        a_d_d_L_(n) = autodiff_d_d_L_(n, a, b, m, t_f)

        @test isapprox(a_d_L_.(t), _d_L_.(t), atol=1e-8, rtol=1e-8)
        @test isapprox(a_d_d_L_.(t), _d_d_L_.(t), atol=1e-8, rtol=1e-8)



        a_d_X = a_d_P_.(t) + a_d_L_.(t)

        @test isapprox(a_d_X, d_X, atol=1e-8, rtol=1e-8)
    end


    # plot!(t,_P_.(t), show=true)
    # plot!(t,_L_.(t), show=true)
    # plot(t,X, show=true)
    # plot!(t,d_X, show=true)
    # plot!(t,d_d_X, show=true)


end