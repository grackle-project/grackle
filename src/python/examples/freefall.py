########################################################################
#
# Free-fall example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import sys
import yt

from pygrackle import \
    chemistry_data, \
    evolve_constant_density, \
    evolve_freefall, \
    setup_fluid_container
from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc
from pygrackle.utilities.testing import \
    dirname
from pygrackle.utilities.model_tests import \
    get_model_set

grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
grackle_data_dir = os.path.join(grackle_install_dir, "input")
output_name = os.path.basename(__file__[:-3]) # strip off ".py"

if __name__=="__main__":
    # If we are running the script through the testing framework,
    # then we will pass in two integers corresponding to the sets
    # of parameters and inputs.
    if len(sys.argv) > 1:
        par_index = int(sys.argv[1])
        input_index = int(sys.argv[2])
        my_chemistry, input_set = get_model_set(
            output_name, par_index, input_index)
        for var, val in input_set.items():
            globals()[var] = val
        output_name = f"{output_name}_{par_index}_{input_index}"

    # Just run the script as is.
    else:
        metallicity = 0.

        # Set solver parameters
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 3
        my_chemistry.metal_cooling = 0
        my_chemistry.dust_chemistry = 0
        my_chemistry.photoelectric_heating = 0
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0
        my_chemistry.CaseBRecombination = 1
        my_chemistry.cie_cooling = 1
        my_chemistry.h2_optical_depth_approximation = 1
        my_chemistry.grackle_data_file = os.path.join(
            grackle_data_dir, "cloudy_metals_2008_3D.h5")

    redshift = 0.

    # Set units
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units  = mass_hydrogen_cgs
    my_chemistry.length_units   = cm_per_mpc
    my_chemistry.time_units     = sec_per_Myr
    my_chemistry.set_velocity_units()

    # set initial density and temperature
    initial_temperature = 50000.
    initial_density     = 1e-1 * mass_hydrogen_cgs # g / cm^3
    final_density       = 1e12 * mass_hydrogen_cgs

    metal_mass_fraction = metallicity * my_chemistry.SolarMetalFractionByMass
    dust_to_gas_ratio = metallicity * my_chemistry.local_dust_to_gas_ratio
    fc = setup_fluid_container(
        my_chemistry,
        density=initial_density,
        temperature=initial_temperature,
        metal_mass_fraction=metal_mass_fraction,
        dust_to_gas_ratio=dust_to_gas_ratio,
        state="ionized",
        converge=False)

    # let the gas cool at constant density from the starting temperature
    # down to a lower temperature to get the species fractions in a
    # reasonable state.
    cooling_temperature = 100.
    data0 = evolve_constant_density(
        fc, final_temperature=cooling_temperature,
        safety_factor=0.1)

    # evolve density and temperature according to free-fall collapse
    data = evolve_freefall(fc, final_density,
                           safety_factor=0.01,
                           include_pressure=True)

    plots = pyplot.loglog(data["density"], data["temperature"],
                          color="black", label="T$_{gas}$")
    if fc.chemistry_data.dust_chemistry == 1:
        plots.extend(
            pyplot.loglog(data["density"], data["dust_temperature"],
                          color="black", linestyle="--", label="T$_{dust}$"))
    pyplot.xlabel("$\\rho$ [g/cm$^{3}$]")
    pyplot.ylabel("T [K]")

    pyplot.twinx()
    plots.extend(
        pyplot.loglog(data["density"], data["H2I_density"] / data["density"],
                      color="red", label="f$_{H2}$"))
    pyplot.ylabel("H$_{2}$ fraction")
    pyplot.legend(plots, [plot.get_label() for plot in plots],
                  loc="lower right")
    pyplot.tight_layout()
    pyplot.savefig(f"{output_name}.png")

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, f"{output_name}.h5", data=data)
