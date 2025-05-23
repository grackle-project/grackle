########################################################################
#
# Cooling rate example script
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
import numpy as np
import os
import sys
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container
from pygrackle.utilities.data_path import grackle_data_dir
from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc
from pygrackle.utilities.model_tests import \
    get_test_variables

output_name = os.path.basename(__file__[:-3]) # strip off ".py"

def gen_plot(fc, data, fname):
    pyplot.loglog(data["temperature"], np.abs(data["cooling_rate"]),
                  color="black")
    pyplot.xlabel('T [K]')
    pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')
    pyplot.tight_layout()
    pyplot.savefig(fname)


if __name__ == "__main__":
    # If we are running the script through the testing framework,
    # then we will pass in two integers corresponding to the sets
    # of parameters and inputs.
    if len(sys.argv) > 1:
        my_vars = get_test_variables(output_name, *sys.argv[1:])
        for var, val in my_vars.items():
            globals()[var] = val

        in_testing_framework = True

    # Just run the script as is.
    else:
        metallicity = 1. # Solar
        redshift = 0.
        specific_heating_rate = 0.
        volumetric_heating_rate = 0.
        # dictionary to store extra information in output dataset
        extra_attrs = {}

        # Set solver parameters
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 0
        my_chemistry.primordial_chemistry = 3
        my_chemistry.metal_cooling = 1
        my_chemistry.UVbackground = 1
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0
        my_chemistry.grackle_data_file = \
          os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

        my_chemistry.use_specific_heating_rate = 1
        my_chemistry.use_volumetric_heating_rate = 1

        in_testing_framework = False

    # max_iterations needs to be increased for the colder temperatures
    my_chemistry.max_iterations = 4*10000

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    my_chemistry.set_velocity_units()

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    metal_mass_fraction = metallicity * my_chemistry.SolarMetalFractionByMass
    temperature = np.logspace(1, 9, 200)
    fc = setup_fluid_container(
        my_chemistry,
        density=mass_hydrogen_cgs,
        temperature=temperature,
        metal_mass_fraction=metal_mass_fraction,
        converge=True)

    if my_chemistry.use_specific_heating_rate:
        fc["specific_heating_rate"][:] = specific_heating_rate
    if my_chemistry.use_volumetric_heating_rate:
        fc["volumetric_heating_rate"][:] = volumetric_heating_rate

    # get data arrays with symbolic units
    data = fc.finalize_data()

    if not in_testing_framework:
        gen_plot(fc, data, fname=f"{output_name}.png")

    yt.save_as_dataset({}, filename=f"{output_name}.h5",
                       data=data, extra_attrs=extra_attrs)
