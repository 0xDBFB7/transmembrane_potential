
##
using Test, ForwardDiff, NumericalIntegration;
using Plots;
using LinearAlgebra;
using BenchmarkTools;
#using Quadmath;

epsilon = 1e-15


# https://m3g.github.io/JuliaNotes.jl/stable/memory/
# https://diffeq.sciml.ai/latest/analysis/sensitivity/

# Nagurka m, Yen V, Benaroya H.
# A Fourier-based method for the suboptimal control of nonlinear dynamical systems.
# Dynamics and Control of Large Structures 1988:77–88.

#const fdcache = ForwardDiffCache()

# https://stackoverflow.com/q/52045626/2179570


# use differentialequations' sensitivity stuff to get the gradient.
#
# Pkg had an issue precompiling diffequationssensitivity. running Pkg.update()
# fixed it.




function L_(t, a, b, m, t_f)
    v = m.*(2*pi / t_f)
    return sum(a.*cos.(v*t)) + sum(b.*sin.(v*t))
end;

# x::Union{Vector,Number}

autodiff_d_L_(t, a, b, m, t_f) = ForwardDiff.derivative(n -> L_(n, a, b, m, t_f), t)
autodiff_d_d_L_(t, a, b, m, t_f) = ForwardDiff.derivative(n -> d_L_(n, a, b, m, t_f), t)
#P
function d_L_(t, a, b, m, t_f)
    v = m.*(2*pi / t_f)
    return sum(v .* -1.0 .* a.*sin.(v*t)) + sum(v .* b.*cos.(v*t))
end;

function d_d_L_(t, a, b, m, t_f)
    v = m.*(2*pi / t_f)
    return sum((v.^2) .* -1.0 .* a.*cos.(v*t)) + sum((v.^2) .* -1.0 .* b.*sin.(2*pi*(m*t)/t_f))
end;

function X_(t,p,a,b,m,t_f)
    return P_(t, p, a, b, m, t_f) + L_(t, a, b, m, t_f)
end

function d_X_(t,p,a,b,m,t_f)
    return d_P_(t, p, a, b, m, t_f) + d_L_(t, a, b, m, t_f)
end

function d_d_X_(t,p,a,b,m,t_f)
    return d_d_P_(t, p, a, b, m, t_f) + d_d_L_(t, a, b, m, t_f)
end

# should the fourier series be always be *zero* at the edges, or simply *the same*?

function X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
    """
    moves polynomial boundary conditions so that the fourier series is
    always zero at the edges for convergence
    """
    return (X_t0 - L_(epsilon, a, b, m, t_f),
            d_X_t0 - d_L_(epsilon, a, b, m, t_f),
            d_d_X_t0 - d_d_L_(epsilon, a, b, m, t_f),

            X_tf - L_(t_f, a, b, m, t_f),
            d_X_tf - d_L_(t_f, a, b, m, t_f),
            d_d_X_tf - d_d_L_(t_f, a, b, m, t_f))
end

function P_BCs_to_p_coefficients(P_BCs, t_f)
    """
    Converts the polynomial boundary conditions that must be enforced into
    the polynomial coefficients
    """
    P_t0, d_P_t0, d_d_P_t0, P_tf, d_P_tf, d_d_P_tf = P_BCs

    p_0 = P_t0
    p_1 = d_P_t0
    p_2 = d_d_P_t0/2
    p_3 = (t_f^2*(-3*d_d_P_t0 + d_d_P_tf) - 4*t_f*(3*d_P_t0 + 2*d_P_tf) - 20*P_t0 + 20*P_tf)/(2*t_f^3)
    p_4 = (t_f^2*(3*d_d_P_t0 - 2*d_d_P_tf)/2 + t_f*(8*d_P_t0 + 7*d_P_tf) + 15*P_t0 - 15*P_tf)/t_f^4
    p_5 = (t_f^2*(-d_d_P_t0 + d_d_P_tf) - 6*t_f*(d_P_t0 + d_P_tf) - 12*P_t0 + 12*P_tf)/(2*t_f^5)

    return (p_0, p_1, p_2, p_3, p_4, p_5)
end

function P_(t, p, a, b, m, t_f)
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = p

    return p_0 + p_1*t + p_2*t^2 + p_3*t^3 + p_4*t^4 + p_5*t^5
end

autodiff_d_P_(t, P_BCs, a, b, m, t_f) = ForwardDiff.derivative(n -> P_(n, P_BCs, a, b, m, t_f), t)
autodiff_d_d_P_(t, P_BCs, a, b, m, t_f) = ForwardDiff.derivative(n -> d_P_(n, P_BCs, a, b, m, t_f), t)

function d_P_(t, p, a, b, m, t_f)
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = p

    return p_1 + 2*p_2*t + 3*p_3*t^2 + 4*p_4*t^3 + 5*p_5*t^4
end

function d_d_P_(t, p, a, b, m, t_f)
    """
    Polynomial from polynomial_system_of_equations.py
    """
    p_0, p_1, p_2, p_3, p_4, p_5 = p

    return 2*p_2 + 6*p_3*t + 12*p_4*t^2 + 20*p_5*t^3
end



"""
Switching back to Python because all of the transmembrane stuff is in that and it might be clunky to
go back and forth with PyCall.

TODO: replace analytic_ as the default; make forwarddiff optional
"""
##