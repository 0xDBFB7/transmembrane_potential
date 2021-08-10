"""

N.B.

Implementations of various polynomial formulae from Yen and Nagurka.
They each appear to be broken in their own ways
because of mistakes in transcription on my part and typos in the papers.

"""



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
