
include("pore_lib.jl")
include("cell_lib.jl")

@testset "coefficient_test" begin
    cell_1 = tl.Cell((0.3), (80), (0.3), (80), (1e-7), (5), (20e-6), (5e-9), py"""np.array([])""")
    l_m_ep = 0.0
    t1 = electroporation_coefficients(cell_1, l_m_ep)
    @test isapprox(t1[1], cell_1.alpha)    
    @test isapprox(t1[2], cell_1.beta)    
    @test isapprox(t1[3], cell_1.gamma)    
    @test isapprox(t1[4], cell_1.phi)    
    @test isapprox(t1[5], cell_1.xi)    

    l_m_ep = 0.29
    cell_2 = tl.Cell((0.3), (80), (0.3), (80), (1e-7 + l_m_ep), (5), (20e-6), (5e-9), py"""np.array([])""")
    t1 = electroporation_coefficients(cell_1, l_m_ep)
    @test isapprox(t1[1], cell_2.alpha, rtol=1e-4, atol=1e-4)    
    @test isapprox(t1[2], cell_2.beta, rtol=1e-4, atol=1e-4)    
    @test isapprox(t1[3], cell_2.gamma, rtol=1e-4, atol=1e-4)    
    @test isapprox(t1[4], cell_2.phi, rtol=1e-4, atol=1e-4)    
    @test isapprox(t1[5], cell_2.xi, rtol=1e-4, atol=1e-4)    

end


# @testset "a" begin

#     t_f = 1e-6
#     a = [2]
#     b = [1.442172e-3]

#     M = 1
#     m = [1.0:M;]
#     t = range(epsilon, stop=t_f, length=500)

#     X_t0 = 0.0
#     X_tf = 1.0
#     d_X_t0 = 0.0
#     d_X_tf = 2.0
#     d_d_X_t0 = 6.1025137
#     d_d_X_tf = -3.4798053


#     P_BCs = X_to_P_BCs(X_t0, d_X_t0, d_d_X_t0, X_tf, d_X_tf, d_d_X_tf, t_f, a, b, m)
#     p = P_BCs_to_p_coefficients(P_BCs, t_f)

#     virus_membrane_thickness = virus.membrane_thickness# 5e-9 # overriding temporarily! TODO: fix the root cause here!
#     cell_v = cell_struct(virus.alpha, virus.beta,virus.gamma,virus.phi,virus.xi,
#                 virus.membrane_permittivity, virus_membrane_thickness, virus.cell_diameter)

#     _X_(n) = X_(n,p,a,b,m,t_f)
#     _d_X_(n) = d_X_(n,p,a,b,m,t_f)
#     _d_d_X_(n) = d_d_X_(n,p,a,b,m,t_f)
#     X = _X_.(t)
#     d_X = _d_X_.(t)
#     d_d_X = _d_d_X_.(t)

#     _d_V_ep(v) = d_V_ep(v, 1e15, cell_v)
#     _d_d_V_ep(v,v2) = d_d_V_ep(v, v2, 1e15, cell_v)

#     @gp t[begin+1:end] diff(_d_V_ep.(X))/(t[end]/500)
#     @gp :- t _d_d_V_ep.(X, d_X)

# end


