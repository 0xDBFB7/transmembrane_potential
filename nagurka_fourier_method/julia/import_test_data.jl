using DelimitedFiles
using Dierckx

function import_talele_test_data()
    pore_density = readdlm("../../test_data/talele/talele_fig_2.csv", ',', skipstart=1, Float64)
    transmembrane_voltage = readdlm("../../test_data/talele/talele_fig_2_2.csv", ',', skipstart=1, Float64)
    
    microsecond = 1e-6
    pore_density_interpolate = Spline1D(pore_density[:,1] * microsecond, pore_density[:,2])
    transmembrane_voltage_interpolate = Spline1D(transmembrane_voltage[:,1] * microsecond, transmembrane_voltage[:,2])
    return pore_density_interpolate, transmembrane_voltage_interpolate
end
