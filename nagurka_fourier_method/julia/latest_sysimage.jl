using PackageCompiler;
@show "d"
using Plots
using Gnuplot
using DifferentialEquations
using Test
@show "d"
using LinearAlgebra
using ForwardDiff
using NumericalIntegration
using BenchmarkTools
@show "d"
using DiffEqSensitivity
using DiffEqParamEstim

pkgs = [:Plots,
        :Gnuplot, 
        :DifferentialEquations,
        :ForwardDiff,
        :NumericalIntegration,
        :BenchmarkTools,
        :DiffEqSensitivity,
        :DiffEqParamEstim
        ]
@show
create_sysimage(pkgs;
                    precompile_execution_file="nagurka_diffeq_lib_test.jl", replace_default=true)

# PackageCompiler.restore_default_sysimage()
