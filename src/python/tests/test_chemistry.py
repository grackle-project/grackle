########################################################################
#
# Chemistry testing functions
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np
import os

from pygrackle import \
    chemistry_data, \
    setup_fluid_container, \
    set_cosmology_units

from pygrackle.utilities.testing import \
    grackle_data_dir, \
    random_logscale, \
    assert_rel_equal, \
    assert_array_less


def test_proper_comoving_units():
    """
    Make sure proper and comoving units systems give the same answer.
    """

    data_file_path = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    my_random_state = np.random.RandomState(7921)
    for current_redshift in [0., 1., 3., 6., 9.]:

        # comoving units
        chem_c = chemistry_data()
        chem_c.use_grackle = 1
        chem_c.with_radiative_cooling = 0
        chem_c.primordial_chemistry = 1
        chem_c.metal_cooling = 1
        chem_c.UVbackground = 1
        chem_c.grackle_data_file = data_file_path
        metal_fraction = 0.1 * chem_c.SolarMetalFractionByMass
        set_cosmology_units(chem_c,
                            current_redshift=current_redshift,
                            initial_redshift=99.)
        fc_c = setup_fluid_container(chem_c, converge=True,
                                     metal_mass_fraction=metal_fraction)
        fc_c.calculate_temperature()
        fc_c.calculate_cooling_time()
        t_sort_c = np.argsort(fc_c["temperature"])
        t_cool_c = fc_c["cooling_time"][t_sort_c] * chem_c.time_units

        # proper units
        chem_p = chemistry_data()
        chem_p.use_grackle = 1
        chem_p.with_radiative_cooling = 0
        chem_p.primordial_chemistry = 1
        chem_p.metal_cooling = 1
        chem_p.UVbackground = 1
        chem_p.grackle_data_file = data_file_path
        chem_p.comoving_coordinates = 0
        chem_p.a_units = 1.0
        chem_p.a_value = 1.0 / (1.0 + current_redshift) / chem_p.a_units
        # Set the proper units to be of similar magnitude to the
        # comoving system to help the solver be more efficient.
        chem_p.density_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.density_units / (1 + current_redshift)**3
        chem_p.length_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.length_units * (1 + current_redshift)
        chem_p.time_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.time_units
        chem_p.velocity_units = chem_p.length_units / chem_p.time_units
        fc_p = setup_fluid_container(chem_p, converge=True,
                                     metal_mass_fraction=metal_fraction)
        fc_p.calculate_temperature()
        fc_p.calculate_cooling_time()
        t_sort_p = np.argsort(fc_p["temperature"])
        t_cool_p = fc_p["cooling_time"][t_sort_p] * chem_p.time_units

        comp = f"\nDU1: {chem_p.density_units:e}, LU1: {chem_p.length_units:e}, " + \
          f"TU1: {chem_p.time_units:e} - DU2: {chem_c.density_units:e}, " + \
          f"LU2: {chem_c.length_units:e}, TU2 {chem_c.time_units:e}."

        rat = t_cool_p / t_cool_c
        err_msg = "Proper and comoving cooling times disagree for " + \
          f"z = {current_redshift} with min/max = {rat.min()}/{rat.max()}." + comp
        assert_rel_equal(t_cool_p, t_cool_c, 4, err_msg=err_msg)


def test_proper_comoving_units_tabular():
    """
    Make sure proper and comoving units systems give the same
    answer with tabular cooling.
    """

    data_file_path = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    my_random_state = np.random.RandomState(19650909)
    for current_redshift in [0., 1., 3., 6., 9.]:

        # comoving units
        chem_c = chemistry_data()
        chem_c.use_grackle = 1
        chem_c.with_radiative_cooling = 0
        chem_c.primordial_chemistry = 0
        chem_c.metal_cooling = 1
        chem_c.UVbackground = 1
        chem_c.grackle_data_file = data_file_path
        set_cosmology_units(chem_c,
                            current_redshift=current_redshift,
                            initial_redshift=99.)
        fc_c = setup_fluid_container(chem_c, converge=False)
        fc_c.calculate_temperature()
        fc_c.calculate_cooling_time()
        t_sort_c = np.argsort(fc_c["temperature"])
        t_cool_c = fc_c["cooling_time"][t_sort_c] * chem_c.time_units

        # proper units
        chem_p = chemistry_data()
        chem_p.use_grackle = 1
        chem_p.with_radiative_cooling = 0
        chem_p.primordial_chemistry = 0
        chem_p.metal_cooling = 1
        chem_p.UVbackground = 1
        chem_p.grackle_data_file = data_file_path
        chem_p.comoving_coordinates = 0
        chem_p.a_units = 1.0
        chem_p.a_value = 1.0 / (1.0 + current_redshift) / chem_p.a_units
        # Set the proper units to be of similar magnitude to the
        # comoving system to help the solver be more efficient.
        chem_p.density_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.density_units / (1 + current_redshift)**3
        chem_p.length_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.length_units * (1 + current_redshift)
        chem_p.time_units = random_logscale(-2, 2, random_state=my_random_state)[0] * \
            chem_c.time_units
        chem_p.velocity_units = chem_p.length_units / chem_p.time_units
        fc_p = setup_fluid_container(chem_p, converge=False)
        fc_p.calculate_temperature()
        fc_p.calculate_cooling_time()
        t_sort_p = np.argsort(fc_p["temperature"])
        t_cool_p = fc_p["cooling_time"][t_sort_p] * chem_p.time_units

        comp = f"\nDU1: {chem_p.density_units:e}, LU1: {chem_p.length_units:e}, " + \
          f"TU1: {chem_p.time_units:e} - DU2: {chem_c.density_units:e}, " + \
          f"LU2: {chem_c.length_units:e}, TU2 {chem_c.time_units:e}."

        rat = t_cool_p / t_cool_c
        err_msg = "Proper and comoving tabules cooling times disagree for " + \
          f"z = {current_redshift} with min/max = {rat.min()}/{rat.max()}." + comp
        assert_rel_equal(t_cool_p, t_cool_c, 4, err_msg=err_msg)


def test_proper_units():
    """
    Make sure two different proper units systems give the same answer.
    """

    data_file_path = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    my_random_state = np.random.RandomState(20150725)
    for current_redshift in [0., 1., 3.]:

        # proper units
        chem_1 = chemistry_data()
        chem_1.use_grackle = 1
        chem_1.with_radiative_cooling = 0
        chem_1.primordial_chemistry = 1
        chem_1.metal_cooling = 1
        chem_1.UVbackground = 1
        chem_1.grackle_data_file = data_file_path
        chem_1.comoving_coordinates = 0
        chem_1.a_units = 1.0
        chem_1.a_value = 1.0 / (1.0 + current_redshift) / chem_1.a_units
        chem_1.density_units = random_logscale(-1, 1, random_state=my_random_state)[0]
        chem_1.length_units = random_logscale(0, 2, random_state=my_random_state)[0]
        chem_1.time_units = random_logscale(0, 2, random_state=my_random_state)[0]
        chem_1.velocity_units = chem_1.length_units / chem_1.time_units
        fc_1 = setup_fluid_container(chem_1, converge=False)
        fc_1.calculate_temperature()
        fc_1.calculate_cooling_time()
        t_sort_1 = np.argsort(fc_1["temperature"])
        t_cool_1 = fc_1["cooling_time"][t_sort_1] * chem_1.time_units

        # proper units
        chem_2 = chemistry_data()
        chem_2.use_grackle = 1
        chem_2.with_radiative_cooling = 0
        chem_2.primordial_chemistry = 1
        chem_2.metal_cooling = 1
        chem_2.UVbackground = 1
        chem_2.grackle_data_file = data_file_path
        chem_2.comoving_coordinates = 0
        chem_2.a_units = 1.0
        chem_2.a_value = 1.0 / (1.0 + current_redshift) / chem_2.a_units
        chem_2.density_units = random_logscale(-28, -26, random_state=my_random_state)[0]
        chem_2.length_units = random_logscale(0, 2, random_state=my_random_state)[0]
        chem_2.time_units = random_logscale(0, 2, random_state=my_random_state)[0]
        chem_2.velocity_units = chem_2.length_units / chem_2.time_units
        fc_2 = setup_fluid_container(chem_2, converge=False)
        fc_2.calculate_temperature()
        fc_2.calculate_cooling_time()
        t_sort_2 = np.argsort(fc_2["temperature"])
        t_cool_2 = fc_2["cooling_time"][t_sort_2] * chem_2.time_units

        comp = f"\nDU1: {chem_1.density_units:e}, LU1: {chem_1.length_units:e}, " + \
          f"TU1: {chem_1.time_units:e} - DU2: {chem_2.density_units:e}, " + \
          f"LU2: {chem_2.length_units:e}, TU2 {chem_2.time_units:e}."

        rat = t_cool_1 / t_cool_2
        err_msg = "Different proper unit system cooling times disagree for " + \
          f"z = {current_redshift} with min/max = {rat.min()}/{rat.max()}." + comp
        assert_rel_equal(t_cool_1, t_cool_2, 4, err_msg=err_msg)


def test_tabulated_mmw_metal_dependence():
    """
    Make sure that increasing the metal mass fraction of a gas increases the
    mean molecular weight, when run in tabulated mode.
    """

    data_file_path = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    my_random_state = np.random.RandomState(723466)
    density_units = random_logscale(-28, -26, random_state=my_random_state)[0]
    length_units = random_logscale(0, 2, random_state=my_random_state)[0]
    time_units = random_logscale(0, 2, random_state=my_random_state)[0]
    velocity_units = length_units / time_units

    current_redshift = 0.

    mmw_vals = []

    for metal_mass_frac in [0.,0.02041]:

        # proper units
        my_chem = chemistry_data()
        my_chem.use_grackle = 1
        my_chem.with_radiative_cooling = 0
        my_chem.primordial_chemistry = 0
        my_chem.metal_cooling = 1
        my_chem.UVbackground = 1
        my_chem.grackle_data_file = data_file_path
        my_chem.comoving_coordinates = 0
        my_chem.a_units = 1.0
        my_chem.a_value = 1.0 / (1.0 + current_redshift) / my_chem.a_units
        my_chem.density_units = density_units
        my_chem.length_units = length_units
        my_chem.time_units = time_units
        my_chem.velocity_units = velocity_units
        fc = setup_fluid_container(my_chem, converge=False,
                                   metal_mass_fraction=metal_mass_frac)
        fc.calculate_mean_molecular_weight()
        mmw_vals.append(fc["mean_molecular_weight"])

    mmw_no_metals, mmw_with_metals = mmw_vals

    assert_array_less(
        mmw_no_metals, mmw_with_metals,
        ("The mmw didn't increase when the metal fraction of the gas was "
         "increased (for primordial_chemisty=0)"))
