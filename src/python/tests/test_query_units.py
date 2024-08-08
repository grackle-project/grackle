########################################################################
#
# Testing the gr_query_units function
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os
import numpy as np
import pytest

from pygrackle import \
    chemistry_data, \
    set_cosmology_units

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from pygrackle.grackle_wrapper import _query_units

_local_dir = os.path.dirname(os.path.abspath(__file__))
def _setup_generic_chemistry_data(initial_redshift, current_redshift = None):
    # construct a generic chemistry_data instance
    # -> it is ONLY set up for comoving coordinates when current_redshift is
    #    not None
    data_file_path = os.sep.join([_local_dir, "..", "..", "..", "input",
                                  "CloudyData_UVB=HM2012.h5"])
    chem = chemistry_data()
    chem.use_grackle = 1
    chem.with_radiative_cooling = 0
    chem.primordial_chemistry = 0
    chem.metal_cooling = 1
    chem.UVbackground = 1
    chem.grackle_data_file = data_file_path
    if current_redshift is not None:
        set_cosmology_units(chem,
                            current_redshift=current_redshift,
                            initial_redshift=initial_redshift)
    else:
        chem.comoving_coordinates = 0
        chem.a_units = 1.0
        chem.a_value = 1.0 / (1.0 + initial_redshift) / chem.a_units
        # Set the proper units to be of similar magnitude to the
        # comoving system to help the solver be more efficient.
        chem.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
        chem.length_units = cm_per_mpc         # 1 Mpc in cm
        chem.time_units = sec_per_Myr          # 1 Myr in s
    chem.initialize()
    return chem


_UNITS_NAMES = ('density_units', 'time_units', 'length_units', 'a_value',
                'a_units', 'velocity_units', 'temperature_units')

def _prefetch_units_vals(chem):
    out = {}
    for name in _UNITS_NAMES:
        if name == 'velocity_units':
            val = chem.get_velocity_units()
        else:
            # note: temperature units is computed on the fly as a property
            val = getattr(chem, name)
        out[name] = val
    return out

@pytest.mark.parametrize("comoving_coordinates,initial_redshift",
                         [(False, 1.0), (True, 3.0), (True, 1.0)])
def test_query_units(comoving_coordinates, initial_redshift):
    # when NOT using comoving coordinates, initial_redshift essentially affects
    # the redshift of the UV background

    # set up a chemistry object
    if comoving_coordinates:
        current_redshift = initial_redshift
    else:
        current_redshift = None
    chem = _setup_generic_chemistry_data(initial_redshift = initial_redshift,
                                         current_redshift = current_redshift)

    # retrieve the initial units-related quantities
    units_at_init = _prefetch_units_vals(chem)

    # first, test retreival of units, when we specify a cosmological scale
    # factor of -1. In this case, we should literally just retrieve the units
    # specified during initialization
    for name in _UNITS_NAMES:
        if _query_units(chem, name, -1) != units_at_init[name]:
            raise AssertionError(f"mismatch when fetching '{name}'")

    
    # next, test retrieval of units, when we specify a cosmological scale
    # factor that exactly matches the initial value. Once again, we should
    # literally just retrieve the units specified during initialization
    a_init = chem.a_value
    for name in _UNITS_NAMES:
        if _query_units(chem, name, a_init) != units_at_init[name]:
            raise AssertionError(f"mismatch when fetching '{name}' with "
                                 "gr_query_units")

    # now, we will test retrieval of units when we specify a cosmological scale
    # factor corresponding to a later redshift
    later_redshift = initial_redshift * 0.5
    a_later = (1.0 / (1.0 + later_redshift)) / chem.a_units
    if not comoving_coordinates:
        # the return value should denote an error if the specified a_value does
        # not EXACTLY match the initial value
        # -> NOTE: this will send messages to stderr
        for name in _UNITS_NAMES:
            assert _query_units(chem, name, a_later) < 0

    else:
        # for the comoving-case, the returned value should match the physical
        # units at the desired cosmological scale_factor
        expected = _prefetch_units_vals(
            _setup_generic_chemistry_data(initial_redshift = initial_redshift,
                                          current_redshift = later_redshift)
        )
        for name in _UNITS_NAMES:
            if name in ('time_units', 'a_units'):
                # these particular quantities should be unchanged
                assert _query_units(chem, name, a_later) == expected[name]
                assert _query_units(chem, name, a_later) == getattr(chem,
                                                                       name)
            else:
                # we don't expect these to be exactly equal
                np.testing.assert_allclose(
                   _query_units(chem, name, a_later), expected[name],
                   rtol = 1e-18, atol = 0,
                   equal_nan = False, # should be no NaNs
                   err_msg = (f"mismatch when fetching '{name}' with "
                              "gr_query_units"))

