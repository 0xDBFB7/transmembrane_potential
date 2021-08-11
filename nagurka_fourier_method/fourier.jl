
using PyCall
push!(pyimport("sys")."path", "../")
 = pyimport("transmembrane_lib")
math.sin(math.pi / 4) # returns ≈ 1/√2 = 0.70710678...


include("nagurka_fourier.jl")
