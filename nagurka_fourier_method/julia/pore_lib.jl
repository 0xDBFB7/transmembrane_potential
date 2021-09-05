
include("nagurka_fourier_lib.jl")
include("cell_lib.jl")


"""

Stop the presses! I've been adding the current term to the capacitor. However, the correct solution is given in Talele:

Once the pores are electrically
induced, the specific membrane conductance increases
because of the increased number of conductive pores, and
can be mathematically modeled as
Iep = Niep = gp Vm (17)
where gp is the specific conductance due to electroporation.
As an effect, the electropermeabilized membrane conductivity sigmap2 is modified to
sigmapp2 =  

where sp is the specific conductivity due to electroporation
and is dynamically altered at each time step as the pore

so instead of just having another term, now the diffeq is bleedin' parametric :(

"""

# scale one equation by the other so to speak? 
# use the chain rule to make the time function dependent on one of the equations? 
# skip the fourier series and go for a polynomial?


function d_pore_density(V_m, N, N0, alpha, q, V_ep)
    # if(abs(V_m) > 2.0)
    #     V_m = 2.0
    # end
    k = (V_m/V_ep)^2
    return alpha * exp(k) * (1.0 - (N/N0)*exp(-q*k))
end;


function electroporation_pore_current(V_m, N, cell)
    """
    The term of the first derivative of the transmembrane potential
    caused by the current flow through the 
    """
    F = 96485.332 # units?
    T = 295.0
    R = 8.314

    r_m = 0.76e-9 # pore radius constant
    pore_solution_conductivity = 13 
    # normally 0.1 mS/cm to S/m, but this peaks out transmembrane at like 8 rather than 3
    # where does 1.3 S/m come from? that's pretty high...
    w0 = 2.65 # ?
    n = 0.15
    v_m = (V_m + epsilon) * (F/(R*T))

    i_ep_term_1 = (pi * (r_m^2) * pore_solution_conductivity * v_m * R * T) / (F * cell.membrane_thickness)

    i_ep_term_2_divisor_1 = exp(v_m)*((w0*exp(w0-n*v_m) - (n*v_m)) / (w0 - (n*v_m)))
    i_ep_term_2_divisor_2 = -((w0*exp(w0+n*v_m) + (n*v_m)) / (w0 + (n*v_m)))
    i_ep_term_2 = exp(v_m - 1) / (i_ep_term_2_divisor_1 + i_ep_term_2_divisor_2)

    i_ep = i_ep_term_1 * i_ep_term_2

    return i_ep * N
end