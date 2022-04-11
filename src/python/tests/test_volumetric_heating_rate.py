########################################################################
#
# Volumetric heating rate tests
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from pygrackle.utilities.testing import \
    assert_almost_equal, \
    random_logscale

def container_setup(density=mass_hydrogen_cgs,
                    current_redshift=0):

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    my_chemistry.primordial_chemistry = 1
    my_chemistry.metal_cooling = 0
    my_chemistry.UVbackground = 0
    my_chemistry.use_volumetric_heating_rate = 1

    # Set units
    my_chemistry.comoving_coordinates = 1
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + current_redshift) / \
      my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs * \
      (1 + current_redshift)**3
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    temperature = np.logspace(1, 6, 51)
    fc = setup_fluid_container(my_chemistry,
                               density=density,
                               temperature=temperature,
                               converge=True)
    fc["volumetric_heating_rate"][:] = 1.e-24 # erg/s/cm^3
    return fc

def get_heating_rate(fc):
    fc.calculate_cooling_time()
    density_proper = fc["density"] / \
        (fc.chemistry_data.a_units *
         fc.chemistry_data.a_value)**(3*fc.chemistry_data.comoving_coordinates)
    cooling_rate = fc.chemistry_data.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper / \
        fc.chemistry_data.a_units**3

    # return the lowest temperature value as it should be all heating
    return cooling_rate[0]

def test_volumetric_heating_rate_density():
    """
    For a constant volumetric heating rate [erg/s/cm^3]
    => net cooling rate [erg*cm^3/s] should go as 1 / rho^2.
    """
    density = mass_hydrogen_cgs
    fc = container_setup(density=density)
    heating_rate = get_heating_rate(fc)

    myrand = np.random.RandomState(seed=8675309)

    for i in range(6):
        new_density = mass_hydrogen_cgs * \
          random_logscale(-1, 1, random_state=myrand)[0]

        for field in fc.density_fields:
            fc[field] *= new_density / density

        new_heating_rate = get_heating_rate(fc)

        assert_almost_equal(new_heating_rate/heating_rate,
                            (density/new_density)**2, decimal=1,
                            err_msg="\nViolates density constraint!")

        density = new_density
        heating_rate = new_heating_rate

def test_volumetric_heating_rate_density_units():
    """
    For a constant volumetric heating rate [erg/s/cm^3]
    => net cooling rate [erg*cm^3/s] should be
       independent of density_units.
    """
    fc = container_setup(density=0.5*mass_hydrogen_cgs)
    density_units = fc.chemistry_data.density_units
    heating_rate = get_heating_rate(fc)

    myrand = np.random.RandomState(seed=8675309)

    for i in range(6):
        new_density_units = density_units * \
          random_logscale(-2, 2, random_state=myrand)[0]
        fc.density_units = new_density_units

        new_heating_rate = get_heating_rate(fc)

        assert_almost_equal(new_heating_rate/heating_rate,
                            1.0, decimal=1,
                            err_msg="\nViolates density_units constraint!")

        density_units = new_density_units
        heating_rate = new_heating_rate

def test_volumetric_heating_rate_redshift():
    """
    For a constant volumetric heating rate [erg/s/cm^3]
    => net cooling rate [erg*cm^3/s] should be
       independent of redshift.
    """
    redshift = 0.
    fc = container_setup(density=0.5*mass_hydrogen_cgs,
                         current_redshift=redshift)
    heating_rate = get_heating_rate(fc)

    myrand = np.random.RandomState(seed=8675309)

    for i in range(6):
        new_redshift = 3 * myrand.random_sample()
        fc.chemistry_data.a_value = 1.0 / (1.0 + new_redshift) / \
          fc.chemistry_data.a_units

        new_heating_rate = get_heating_rate(fc)

        assert_almost_equal(new_heating_rate/heating_rate,
                            1.0, decimal=1,
                            err_msg="\nViolates redshift constraint!")

        redshift = new_redshift
        heating_rate = new_heating_rate

def test_volumetric_heating_rate_a_units():
    """
    For a constant volumetric heating rate [erg/s/cm^3]
    => net cooling rate [erg*cm^3/s] should be
       independent of a_units.
    """
    redshift = 0.
    fc = container_setup(density=0.5*mass_hydrogen_cgs,
                         current_redshift=redshift)
    heating_rate = get_heating_rate(fc)
    velocity_units = fc.chemistry_data.velocity_units

    myrand = np.random.RandomState(seed=8675309)

    for i in range(6):
        fc.chemistry_data.a_units = \
          random_logscale(-0.5, 0.5, random_state=myrand)[0]
        fc.chemistry_data.a_value = 1.0 / (1.0 + redshift) / \
          fc.chemistry_data.a_units
        new_velocity_units = fc.chemistry_data.a_units * \
          (fc.chemistry_data.length_units / fc.chemistry_data.a_value) / \
          fc.chemistry_data.time_units
        fc.chemistry_data.velocity_units = new_velocity_units
        fc["energy"] *= (velocity_units / new_velocity_units)**2

        new_heating_rate = get_heating_rate(fc)

        assert_almost_equal(new_heating_rate/heating_rate,
                            1.0, decimal=1,
                            err_msg="\nViolates a_units constraint!")

        velocity_units = new_velocity_units
        heating_rate = new_heating_rate
