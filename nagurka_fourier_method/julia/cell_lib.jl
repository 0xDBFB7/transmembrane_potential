using PyCall
pushfirst!(PyVector(pyimport("sys")["path"]), "../../")
tl = pyimport("transmembrane_lib")

py"""
import numpy as np
"""

virus = tl.default_virus(py"""np.array([])""")
host_cell = tl.default_host_cell(py"""np.array([])""")

# this darn redundancy needed because virus is global pystruct and slow mutable
struct cell_struct 
    alpha
    beta
    gamma
    phi
    xi
    membrane_permittivity
    membrane_thickness
    cell_diameter
end
