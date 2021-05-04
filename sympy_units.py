from sympy import symbols
from sympy.physics.units.systems import SI
from sympy.physics.units import length, mass, acceleration, force
from sympy.physics.units import gravitational_constant as G
from sympy.physics.units.systems.si import dimsys_SI
import sympy.physics.units as units
#import pprint as pps
import sympy
import sympy.physics.units.util as util
from dataclasses import dataclass
from sympy.simplify.radsimp import collect
from sympy.assumptions.refine import refine
from sympy import init_printing
from sympy.simplify.powsimp import powsimp
init_printing()
'''

To get a better idea of the possible nondimensionalization,
want to get the units of the 6 coefficients.

'''

# pp = pps.PrettyPrinter(indent=4)

t0 = 0
tstop = 0.05e-6
dt = 0.00005e-9

eps0 = units.farads / units.meter
m = units.meters

Spm = units.siemens / units.meter




#changes needed for sympy
@dataclass
class Cell:
    extracellular_conductivity: 'typing.Any' # S/m
    extracellular_permittivity: 'typing.Any' # relative
    intracellular_conductivity: 'typing.Any' # S/m
    intracellular_permittivity: 'typing.Any' # relative
    membrane_conductivity: 'typing.Any' # S/m
    membrane_permittivity: 'typing.Any' # relative
    cell_diameter: 'typing.Any' # meters
    membrane_thickness: 'typing.Any'


    def __post_init__(self):
        e_o, e_i, e_m, R, l_o, l_i, l_m, d, a_1, a_2, a_3, b_1, b_2, b_3 = symbols('e_o e_i e_m R l_o l_i l_m d a_1 a_2 a_3 b_1 b_2 b_3')
        e_o = self.extracellular_permittivity * eps0 # S/m
        e_i = self.intracellular_permittivity * eps0 #S/m
        e_m = self.membrane_permittivity * eps0 #S/m
        R = self.cell_diameter / 2
        self.R = R

        l_o = self.extracellular_conductivity*Spm # S/m
        l_i = self.intracellular_conductivity*Spm #S/m
        l_m = self.membrane_conductivity*Spm #S/m

        d = self.membrane_thickness
        # epsilon_0

        sub1, sub2 = symbols('sub1 sub2')
        sub1 = (3 * (R**2) - 3 * d * R + d**2)
        sub2 = (3 * d * R - d**2)

        a_1 = 3 * d * l_o * ((l_i * (sub1)) + l_m*(sub2)) #eq.9a
        a_2 = 3 * d * ((l_i * e_o + l_o * e_i) * sub1 + (l_m * e_o + l_o * e_m) * sub2)
        a_3 = 3 * d * e_o * (e_i * (sub1) + e_m * sub2)

        b_1 = 2 * R**3 * (l_m +     2*l_o) * (l_m + (1/2) * l_i) + 2 * (R-d)**3 * (l_m - l_o) * (l_i - l_m)

        b_2 = 2 * R**3 * (l_i * ((1/2) * e_m + e_o) + l_m * ((1/2)*e_i + 2*e_m + 2*e_o) + l_o * (e_i + 2 * e_m)) + (2 * (R - d)**3\
        * (l_i * (e_m - e_o) + l_m * (e_i - 2*e_m + e_o) - l_o * (e_i - e_m))) # is this truly a multiply, or a cross?


        b_3 = 2 * R**3 * (e_m + 2*e_o) * (e_m + (1/2) * e_i) + 2 * (R-d)**3 * (e_m - e_o) * (e_i - e_m)

        sympy.pprint(powsimp(util.quantity_simplify(refine(a_1.simplify()))).as_coeff_Mul()[1])
        print()
        sympy.pprint(powsimp(util.quantity_simplify(refine(a_2.simplify()))).as_coeff_Mul()[1])
        print()
        sympy.pprint(powsimp(util.quantity_simplify(refine(a_3.simplify()))).as_coeff_Mul()[1])
        print()
        sympy.pprint(powsimp(util.quantity_simplify(refine(b_1.simplify()))).as_coeff_Mul()[1])
        print()
        sympy.pprint(powsimp(util.quantity_simplify(refine(b_2.simplify()))).as_coeff_Mul()[1])
        print()
        sympy.pprint(powsimp(util.quantity_simplify(refine(b_3.simplify()))).as_coeff_Mul()[1])
        print()
        print()
        # print()
        # sympy.pprint(dimsys_SI.get_dimensional_dependencies(powsimp(util.quantity_simplify(refine(a_1.simplify()))).as_coeff_Mul()[1]))
        # print()

host_cell = Cell(0.3, 80, 0.3, 80, 1e-7, 5, 1*m, 1*m)

# virus = Cell(0.3, 80, 0.005, 30, 1e-8, 80, 50e-9*m, 14e-9*m)

# print(pp.pprint(virus.__dict__))
