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
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

if __name__ == "__main__":
    current_redshift = 0.

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    if 'PRIMORDIAL_CHEM' in os.environ:
        my_chemistry.primordial_chemistry = int(os.environ['PRIMORDIAL_CHEM'])
    else:
        my_chemistry.primordial_chemistry = 3
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    my_chemistry.grackle_data_file = os.sep.join(
        [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"])

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    temperature = np.logspace(1, 9, 200)
    fc = setup_fluid_container(my_chemistry,
                               temperature=temperature,
                               converge=True)

    fc.calculate_temperature()
    fc.calculate_cooling_time()

    density_proper = fc["density"] / \
        (my_chemistry.a_units *
         my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)
    cooling_rate = fc.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper

    data = {}
    t_sort = np.argsort(fc["temperature"])
    for field in fc.density_fields:
        data[field] = yt.YTArray(fc[field][t_sort] *
                                 my_chemistry.density_units, "g/cm**3")
    data["energy"]       = yt.YTArray(fc["energy"][t_sort], "erg/g")
    data["temperature"]  = yt.YTArray(fc["temperature"][t_sort], "K")
    data["pressure"]     = yt.YTArray(fc["pressure"][t_sort], "dyne/cm**2")
    data["cooling_time"] = yt.YTArray(fc["cooling_time"][t_sort], "s")
    data["cooling_rate"] = yt.YTArray(cooling_rate[t_sort], "erg*cm**3/s")

    pyplot.loglog(data["temperature"], data["cooling_rate"],
                  color="black")
    pyplot.xlabel('T [K]')
    pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')

    # save data arrays as a yt dataset
    if 'PRIMORDIAL_CHEM' in os.environ:
        ds_name = 'cooling_rate.pc%s.h5' % os.environ['PRIMORDIAL_CHEM']
        im_name = 'cooling_rate.pc%s.png' % os.environ['PRIMORDIAL_CHEM']
    else:
        ds_name = 'cooling_rate.h5'
        im_name = 'cooling_rate.png'
    pyplot.savefig(im_name)
    yt.save_as_dataset({}, ds_name, data)
