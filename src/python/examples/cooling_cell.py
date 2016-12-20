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
import yt

from pygrackle import \
    FluidContainer, \
    chemistry_data, \
    evolve_constant_density

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

tiny_number = 1e-20

if __name__ == "__main__":
    current_redshift = 0.

    # Set initial values
    density             = 0.1 # g /cm^3
    initial_temperature = 1.e6 # K
    final_time          = 100.0 # Myr

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 0
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    grackle_dir = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    my_chemistry.grackle_data_file = os.sep.join(
        [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"])

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Myr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = density
    if my_chemistry.primordial_chemistry > 0:
        fc["HI"][:] = 0.76 * fc["density"]
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - 0.76) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 1:
        fc["H2I"][:] = tiny_number * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
        fc["de"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = 0.1 * fc["density"] * \
          my_chemistry.SolarMetalFractionByMass

    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc["energy"][:] = initial_temperature / \
        fc.chemistry_data.temperature_units
    fc.calculate_temperature()
    fc["energy"][:] *= initial_temperature / fc["temperature"]

    # timestepping safety factor
    safety_factor = 0.01

    # let gas cool at constant density
    data = evolve_constant_density(
        fc, final_time=final_time,
        safety_factor=safety_factor)

    p1, = pyplot.loglog(data["time"].to("Myr"), data["temperature"],
                        color="black", label="T")
    pyplot.xlabel("Time [Myr]")
    pyplot.ylabel("T [K]")

    data["mu"] = data["temperature"] / \
        (data["energy"] * (my_chemistry.Gamma - 1.) *
         fc.chemistry_data.temperature_units)
    pyplot.twinx()
    p2, = pyplot.semilogx(data["time"].to("Myr"), data["mu"],
                          color="red", label="$\\mu$")
    pyplot.ylabel("$\\mu$")
    pyplot.legend([p1,p2],["T","$\\mu$"], fancybox=True,
                  loc="center left")
    pyplot.savefig("cooling_cell.png")

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, "cooling_cell.h5", data)
