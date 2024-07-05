########################################################################
#
# Cooling cell example script
#
#  This will initialize a single cell at a given temperature,
#  iterate the cooling solver for a fixed time, and output the
#  temperature vs. time.
#
#
# Copyright (c) 2015-2016, Grackle Development Team.
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

if __name__ == "__main__":
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
        metallicity = 0.1 # Solar
        redshift = 0.

        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 0
        my_chemistry.metal_cooling = 1
        my_chemistry.UVbackground = 1
        my_chemistry.grackle_data_file = \
          os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    density = 0.1 * mass_hydrogen_cgs # g /cm^3
    temperature = 1e6 # K
    final_time = 100. # Myr

    # Set units
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs
    my_chemistry.length_units = cm_per_mpc
    my_chemistry.time_units = sec_per_Myr
    my_chemistry.set_velocity_units()

    metal_mass_fraction = metallicity * my_chemistry.SolarMetalFractionByMass
    fc = setup_fluid_container(
        my_chemistry,
        density=density,
        temperature=temperature,
        metal_mass_fraction=metal_mass_fraction,
        state="ionized",
        converge=True)

    # evolve gas at constant density
    data = evolve_constant_density(
        fc, final_time=final_time,
        safety_factor=0.01)

    p1, = pyplot.loglog(data["time"].to("Myr"),
                        data["temperature"],
                        color="black", label="T")
    pyplot.xlabel("Time [Myr]")
    pyplot.ylabel("T [K]")
    pyplot.twinx()
    p2, = pyplot.semilogx(data["time"].to("Myr"),
                          data["mean_molecular_weight"],
                          color="red", label="$\\mu$")
    pyplot.ylabel("$\\mu$")
    pyplot.legend([p1,p2],["T","$\\mu$"], fancybox=True,
                  loc="center left")
    pyplot.tight_layout()
    pyplot.savefig(f"{output_name}.png")

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, f"{output_name}.h5", data=data)
