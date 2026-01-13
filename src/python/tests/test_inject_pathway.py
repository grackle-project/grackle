########################################################################
#
# Basic tests that the python layer has appropriate access to relevant
# injection pathway information. Checks of the correctness of the
# injection pathway information are handled in the C++ test-suite.
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os.path

from gracklepy import chemistry_data, FluidContainer
from gracklepy.utilities.data_path import grackle_data_dir
from gracklepy.utilities.physical_constants import (
    mass_hydrogen_cgs,
    sec_per_Myr,
    cm_per_mpc,
)


def test_info_is_accessible():
    # the main objective here is to make sure that the relevant information
    # can be accessed by the python layer
    # -> the primary concern here is that we might change a convention about
    #    *HOW* the information is accessed and forget to propogate the change
    #    into the python layer

    # first, lets set up a chemistry_data instance with settings that should
    # *always* involve multiple grain-species fields
    # - the details here will slightly change with time
    # - we may want to look into configuring presets like the C++ tests

    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 4
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.grackle_data_file = os.path.join(
        grackle_data_dir, "CloudyData_UVB=HM2012.h5"
    )

    my_chemistry.dust_chemistry = 1
    my_chemistry.metal_chemistry = 1
    my_chemistry.dust_species = 1
    my_chemistry.use_dust_density_field = 1
    my_chemistry.multi_metals = 1

    # Set units
    redshift = 0.0
    my_chemistry.comoving_coordinates = 0
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + redshift) / my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs
    my_chemistry.length_units = cm_per_mpc
    my_chemistry.time_units = sec_per_Myr
    my_chemistry.set_velocity_units()

    # initialize my_chemistry
    assert my_chemistry.initialize() == 1

    # let's confirm that we properly require injection pathway fields
    fc = FluidContainer(my_chemistry, n_vals=1)
    inject_pathway_yield_density_names = fc.inject_pathway_density_yield_fields
    n_pathways = len(inject_pathway_yield_density_names)
    assert n_pathways > 0, (
        "internal logic appears to have broken for accessing the modelled injection "
        "pathway names from the python layer"
    )

    # let's confirm that we properly access gas nuclide yields
    nuclide_gas_yields = my_chemistry._experimental_nuclide_gas_inj_path_yields()
    assert len(nuclide_gas_yields) > 0, (
        "internal logic appears to have broken for accessing (from the gracklepy "
        "layer) the assumed gas-phase yields of each relevant metal nuclide"
    )
    assert "C" in nuclide_gas_yields  # extra sanity check!

    # let's confirm that we properly list grain nuclide yields
    grain_yields = my_chemistry._experimental_grain_inj_path_yields()
    assert len(grain_yields) > 0, (
        "internal logic appears to have broken for accessing (from the gracklepy "
        "layer) the assumed yields of each relevant grain-species"
    )
    assert "AC_dust" in grain_yields  # extra sanity check!
