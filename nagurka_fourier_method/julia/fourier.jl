include("nagurka_fourier_lib.jl")

using PyCall
push!(pyimport("sys")."path", "../../")
tl = pyimport("transmembrane_lib")

using Gnuplot

pore_N0, pore_alpha, pore_q, pore_V_ep = tl.pore_N0, tl.pore_alpha, tl.pore_q, tl.pore_V_ep

# https://diffeqparamestim.sciml.ai/dev/methods/optimization_based_methods/
# might be useful

# or, more directly, https://diffeqparamestim.sciml.ai/dev/

# https://gcalderone.github.io/Gnuplot.jl/v1.1.0/

# http://julianlsolvers.github.io/Optim.jl/v0.9.3/algo/nelder_mead/


function 
k = (V_m/V_ep)**2.0

return alpha * np.exp(k) * (1.0 - (N/N0)*np.exp(-q*k))

M = 3

function f(dX,X,p,t)
    a = p[1:M]
    b = p[1:M]

    dX[1] = p[1]*u[1] - u[1]*u[2]
    dX[2] = -3*u[2] + u[1]*u[2]
end

@show pore_N0

#@gp t "with lines"
