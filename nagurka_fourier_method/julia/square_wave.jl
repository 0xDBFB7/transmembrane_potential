
"""
Switching back to Python because all of the transmembrane stuff is in that and it might be clunky to
go back and forth with PyCall.

TODO: replace analytic_ as the default; make forwarddiff optional
"""
##

function square_wave(M)
    # M=30
    a = (zeros(M))
    b = (zeros(M))
    for i in [1.0:2.0:M;]
        b[Int(i)] = (1/(i))
    end
     return a, b 
 end 

 
function init_step_function(params, E, edge_rise_time, start_time)

    edge_rise_time = 1e-9
    k = 1 / (edge_rise_time / params.T0)
    x0 = start_time / params.T0

    peak = E * (params.cell_h.gamma / params.cell_h.xi) 
    # really should do one cell at a time. 
    # Most users of this code won't want to do exactly two cells.

    control_function(t) = evaluate_step_function(t, peak, k, x0)
    return control_function
end

function evaluate_step_function(t, peak, k, x0)
    v = logistic_curve( t, peak, k, x0)
    d_v = d_logistic_curve( t, peak, k, x0)
    d_d_v = d_d_logistic_curve( t, peak, k, x0)
    return v, d_v, d_d_v
end


function logistic_curve(x, peak, k, k0)
    # Nice smooth step function.
    # https://calculus.subwiki.org/wiki/Logistic_function
    return peak / (1+exp(-(k * (x-k0))))
end

function d_logistic_curve(x, peak, k, k0)
    # return logistic_curve(x, peak, k, k0) * (peak - logistic_curve(x, peak, k, k0))
    return k.*peak.*exp(-k.*(-k0 + x))./(1 + exp(-k.*(-k0 + x))).^2
end

function d_d_logistic_curve(x, peak, k, k0)
    return -k.^2 .* peak.*(1 - 2*exp(k.*(k0 - x))./(exp(k.*(k0 - x)) + 1)).*exp(k.*(k0 - x))./(exp(k.*(k0 - x)) + 1).^2
    # return logistic_curve(x, peak, k, k0) * (peak - logistic_curve(x, peak, k, k0)) * (peak - 2*logistic_curve(x, peak, k, k0))
end
