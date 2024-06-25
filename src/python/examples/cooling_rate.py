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
    my_dir = os.path.dirname(os.path.abspath(__file__))
    grackle_data_file = os.path.join(
        my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012.h5")
    my_chemistry.grackle_data_file = grackle_data_file

    my_chemistry.use_specific_heating_rate = 1
    my_chemistry.use_volumetric_heating_rate = 1

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    my_chemistry.set_velocity_units()

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    temperature = np.logspace(1, 9, 200)
    fc = setup_fluid_container(
        my_chemistry,
        temperature=temperature,
        metal_mass_fraction=my_chemistry.SolarMetalFractionByMass,
        converge=True)

    fc["specific_heating_rate"][:] = 0.
    fc["volumetric_heating_rate"][:] = 0.

    # get data arrays with symbolic units
    data = fc.finalize_data()

    pyplot.loglog(data["temperature"], np.abs(data["cooling_rate"]),
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
    pyplot.tight_layout()
    pyplot.savefig(im_name)
    yt.save_as_dataset({}, ds_name, data)
