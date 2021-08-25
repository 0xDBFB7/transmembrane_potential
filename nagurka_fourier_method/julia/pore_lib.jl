
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

"""

function d_pore_density(V_m, N, N0, alpha, q, V_ep)
    # if(abs(V_m) > 2.0)
    #     V_m = 2.0
    # end
    k = (V_m/V_ep)^2
    return alpha * exp(k) * (1.0 - (N/N0)*exp(-q*k))
end;


function d_V_ep(V_m, N, cell)
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
    # I = C_m dv/dt
    # dv_dt = I/C_m
    # permittivity already has the eps0 in it 
    A = 4*pi*(cell.cell_diameter/2)^2
    C_m = cell.membrane_permittivity * A / cell.membrane_thickness 
    d_V_ep = -(i_ep * N / C_m)

    # C_m and A should maybe go in transmembrane_lib

    return d_V_ep

end

function d_d_V_ep(V_m, d_V_m, N, cell)
    """
    Derived via nagurka_math_2.py, then sagemath Maxima backend
    Here the V_m -> v_m normalization is baked in!
    """
    
    F = 96485.332 # units?
    T = 295.0
    R = 8.314

    r_m = 0.76e-9 # pore radius constant
    sigma = pore_solution_conductivity = 13 
    # normally 0.1 mS/cm to S/m, but this peaks out transmembrane at like 8 rather than 3
    # where does 1.3 S/m come from? that's pretty high...
    w0 = 2.65 # ?
    n = 0.15
    diameter = cell.cell_diameter
    e = Float64(MathConstants.e)
    k = cell.membrane_permittivity

    v_m = (V_m + epsilon) * (F/(R*T))
    d_v_m = (d_V_m + epsilon) * (F/(R*T))

    d_d_V_ep = (N*R*T*r_m^2*sigma*e^(v_m - 1)*v_m*d_v_m/(F*diameter^2*k*((w0*e^(-n*v_m + w0) - n*v_m)*e^v_m/(n*v_m - w0) 
                + (w0*e^(n*v_m + w0) + n*v_m)/(n*v_m + w0))) - ((w0*e^(-n*v_m + w0) - n*v_m)*e^v_m*d_v_m/(n*v_m - w0) 
                - (w0*e^(-n*v_m + w0) - n*v_m)*n*e^v_m*d_v_m/(n*v_m - w0)^2 - (n*w0*e^(-n*v_m + w0)*d_v_m 
                + n*d_v_m)*e^v_m/(n*v_m - w0) - (w0*e^(n*v_m + w0) + n*v_m)*n*d_v_m/(n*v_m + w0)^2 + (n*w0*e^(n*v_m 
                + w0)*d_v_m + n*d_v_m)/(n*v_m + w0))*N*R*T*r_m^2*sigma*e^(v_m - 1)*v_m/(F*diameter^2*k*((w0*e^(-n*v_m + w0) 
                - n*v_m)*e^v_m/(n*v_m - w0) + (w0*e^(n*v_m + w0) + n*v_m)/(n*v_m + w0))^2) + N*R*T*r_m^2*sigma*e^(v_m 
                - 1)*d_v_m/(F*diameter^2*k*((w0*e^(-n*v_m + w0) - n*v_m)*e^v_m/(n*v_m - w0) 
                + (w0*e^(n*v_m + w0) + n*v_m)/(n*v_m + w0))))

    return d_d_V_ep

end