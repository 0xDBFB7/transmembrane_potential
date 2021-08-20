include("nagurka_diffeq_lib.jl")

M = 3
m = [1.0:M;]
t_f = 1e-6

a = [0.0, 0.0, 1.0]
b = [0.0, 0.0, 0.1]

c = zeros(6)

O = zeros((2*M)+6)
params = transmembrane_params(O, t_f, M, m)
params.O[1:M] = a
params.O[M+1:(2*M)] = b
params.O[(2*M)+1:(2*M)+6] = c

d = zeros(length(instances(svars)))
s = zeros(length(instances(svars)))
t=epsilon 
@profview [transmembrane_diffeq(d,s,params,i) for i in range(0.0,stop=1.0,length=100)]
# not super helpful for this use