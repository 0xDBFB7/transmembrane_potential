using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "../../")
tl = pyimport("transmembrane_lib")

py"""
import numpy as np
"""



virus = tl.default_virus(py"""np.array([])""")
host_cell = tl.default_host_cell(py"""np.array([])""")

function py_cell_to_julia_struct(input_cell)
    return cell_struct(
        input_cell.alpha,
        input_cell.beta,
        input_cell.gamma,
        input_cell.phi,
        input_cell.xi,
        input_cell.alpha_ep,
        input_cell.beta_ep,
        input_cell.gamma_ep,
        input_cell.phi_ep,
        input_cell.xi_ep,
        input_cell.membrane_permittivity,
        input_cell.membrane_thickness,
        input_cell.pore_solution_conductivity,
        input_cell.cell_diameter,
        input_cell.R)
end

# this darn redundancy needed because virus is global pystruct and slow mutable
struct cell_struct 
    alpha
    beta
    gamma
    phi
    xi
    alpha_ep
    beta_ep
    gamma_ep
    phi_ep
    xi_ep
    membrane_permittivity
    membrane_thickness
    pore_solution_conductivity
    cell_diameter
    R
end
