########################################################################
#
# testing hdf5 dumps
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os
import tempfile
import numpy as np
import pytest

from pygrackle import FluidContainer
from pygrackle.debug_tools import load_FluidContainer_from_h5dump
from pygrackle.grackle_wrapper import _h5dump_state
from pygrackle.utilities.init_from_dicts import \
    _configure_chemistry_data, \
    _configure_codeunits
from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc
from pygrackle.utilities.testing import assert_allequal_arraydict
from pygrackle.utilities.test_scenarios import grid_initial_values
from pygrackle.utilities.units import set_cosmology_units


def _initial_props(primordial_chemistry = 0, initial_redshift = 0.5,
                   comoving = False, use_initial_units = True):
    """
    This function sorta sets up the testing scenario.

    This always, sets up and initializes the chemistry_data instance
    and a dict of fields to be used in a test problem.

    In the case of a cosmological simulation, there is the additional option
    to change the fields.
    """
    # determine initial "code_units" & the code_units used during the calc
    if comoving:
        initial_code_units = set_cosmology_units(
            None, current_redshift = initial_redshift,
            initial_redshift = initial_redshift)
        if use_initial_units:
            cur_code_units = initial_code_units.copy()
        else:
            cur_code_units = set_cosmology_units(
                None, current_redshift = 0.0,
                initial_redshift = initial_redshift)
    else:
        a_units = 1.0
        initial_code_units = dict(
            comoving_coordinates = 0, # proper units
            a_units = a_units,
            a_value = 1. / (1. + initial_redshift) / a_units,
            density_units = mass_hydrogen_cgs, # rho = 1.0 is 1.67e-24 g
            length_units = cm_per_mpc,         # 1 Mpc in cm
            time_units = sec_per_Myr,          # 1 Myr in s
        )
        if use_initial_units:
            cur_code_units = initial_code_units.copy()
        else:
            raise RuntimeError("Must use initial units in non-cosmo case")

    # determine remaining chemistry_parameters
    my_dir = os.path.dirname(os.path.abspath(__file__))
    metal_cooling = 1
    chemistry_parameters = dict(
        use_grackle = 1,
        with_radiative_cooling = 1,
        primordial_chemistry = primordial_chemistry,
        metal_cooling = metal_cooling,
        UVbackground = 1,
        self_shielding_method = 0,
        H2_self_shielding = 0,
        grackle_data_file = os.path.join(
            my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012.h5")
    )

    # now extract the field data used in the calculation
    initial_fluid_config = grid_initial_values(
        primordial_chem_flag = primordial_chemistry,
        metal_cooling_flag = metal_cooling,
        # make sure to initialize with the code-units that are actually used
        # during the calculation (if we used the initialization-code-units, we
        # would need to convert the units on all of the fields)
        code_units = cur_code_units,
        target_T_cgs_arr = [1e4, 1e5, 1e6, 1e7],
        target_Hmassdensity_cgs_arr = mass_hydrogen_cgs * np.array([0.01, 100]),
        # use default metallicity
    )

    # currently, we can only handle 1D arrays...
    for k in initial_fluid_config:
        initial_fluid_config[k].shape = (-1,)

    return (chemistry_parameters, initial_code_units, cur_code_units,
            initial_fluid_config)

def _test_h5dump(h5dump_fname, primordial_chemistry = 0, comoving = False):
    # essentially, we want to test a round-trip
    # -> we setup fully setup Grackle for a test-problem
    # -> before any calculations we dump all of the Grackle configuration and
    #    the test problem from disk
    # -> next, we carry out the actual test-calculation with the original
    #    Grackle-configuration
    # -> then, we load in the dumped data, reconfigure Grackle and setup the
    #    field data based on just what is in the dump and run the
    #    test-calculation with that setup
    # -> finally, we compare the results from both test-calculations. The
    #    results should be bitwise identical.

    tmp = _initial_props(primordial_chemistry = primordial_chemistry,
                         initial_redshift = 0.5,
                         comoving = comoving,
                         use_initial_units = True)
    chemistry_parameters, initial_code_u, cur_code_u = tmp[:-1]
    initial_fluid_config = tmp[-1]

    dt = 1000 * 3.156e7 / cur_code_u['time_units'] # <- 1 kyr

    results_l = []

    # run the test-calculation twice
    for i in range(2):
        # first, fully initialize grackle, update units, record initial
        #        field-data, and prepare FluidContainer
        if i == 0: # this is the first time we set things up
            my_chem = _configure_chemistry_data(chemistry_parameters,
                                                initial_code_u)
            _configure_codeunits(my_chem, cur_code_u) # overwrite code units

            fc = FluidContainer(my_chem,
                                n_vals = initial_fluid_config["density"].size)
            for k,v in initial_fluid_config.items():
                fc[k][...] = v

            # before we go any further, let's dump to hdf5
            initial_code_u["velocity_units"] = (initial_code_u["length_units"] /
                                                initial_code_u["time_units"])
            _h5dump_state(fc, fname = h5dump_fname,
                          initial_code_units = initial_code_u)
        else: # here, (the second time) - we load the dumped data
            fc = load_FluidContainer_from_h5dump(h5dump_fname)

        # now, solve chemistry
        fc.solve_chemistry(dt)
        results_l.append(fc)

    assert_allequal_arraydict(
        actual = results_l[1], reference = results_l[0],
        err_msg = (
            'Comparison failure for h5dump-debugger utility. In this scenario, '
            'we dumped all grackle data to disk before ANY calculations. Next, '
            'we solved chemistry with the initial configuration (the reference '
            'answers). Then, we used the data dump to reinitialize Grackle '
            'from scratch and used that setup to solve chemistry. The results '
            'should be bitwise identical.'))

@pytest.mark.parametrize("comoving", [False, True])
@pytest.mark.parametrize("primordial_chemistry", [0,1,2,3])
def test_h5dump(primordial_chemistry, comoving):
    with tempfile.TemporaryDirectory() as tmpdirname:
        path = os.path.join(tmpdirname, 'my-test-dump.h5')
        _test_h5dump(path, primordial_chemistry = primordial_chemistry,
                     comoving = comoving)
