"""
N.B.

Attempt to use Symbolics.jl to obtain the solution to the differential equation.
And a test of PackageCompiler to speed up JIT compilation.
This is just the Symbolics.jl tutorial file.

"""


using ModelingToolkit, OrdinaryDiffEq
using Plots;
gr(show = true)


# Unfortunately Symbolics.jl doesn't seem to support symbolic ODE solving at this point.

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(D(x)) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

sys = ODESystem(eqs)
sys = ode_order_lowering(sys)

print(sys)

u0 = [D(x) => 2.0,
      x => 1.0,
      y => 0.0,
      z => 0.0]

p  = [σ => 28.0,
      ρ => 10.0,
      β => 8/3]

# sum(p->p.Q_flow, ps)

tspan = (0.0,100.0)
prob = ODEProblem(sys,u0)
sol = solve(prob,Tsit5())
plot(sol,vars=(x,y))
