########################################################################
#
# yt fields for grackle functions
#
#
# Copyright (c) Enzo/Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os
import sys
import yt

from pygrackle import add_grackle_fields
from pygrackle.utilities.data_path import grackle_data_dir
from pygrackle.utilities.model_tests import \
    get_test_variables

output_name = os.path.basename(__file__[:-3]) # strip off ".py"

DS_NAME = "IsolatedGalaxy/galaxy0030/galaxy0030"
ds_path = os.path.join(os.environ.get('YT_DATA_DIR', default='.'), DS_NAME)

if __name__ == "__main__":
    # If we are running the script through the testing framework,
    # then we will pass in two integers corresponding to the sets
    # of parameters and inputs.
    if len(sys.argv) > 1:
        my_vars = get_test_variables(output_name, *sys.argv[1:])
        for var, val in my_vars.items():
            globals()[var] = val

        if 'YT_DATA_DIR' not in os.environ:
            raise RuntimeError(
                "YT_DATA_DIR env var must be defined when called by test suite"
            )

    # Just run the script as is.
    else:
        # dictionary to store extra information in output dataset
        extra_attrs = {}

    ds = yt.load(ds_path)

    grackle_data_file = os.path.join(grackle_data_dir, "CloudyData_UVB=HM2012.h5")

    grackle_pars = {'grackle_data_file': grackle_data_file,
                    'UVbackground': 1,
                    'h2_on_dust': 1}

    add_grackle_fields(ds, parameters=grackle_pars)

    sp = ds.sphere(ds.domain_center, (10, 'kpc'))

    fields = [
        ("gas", "grackle_cooling_time"),
        ("gas", "grackle_gamma"),
        ("gas", "grackle_mean_molecular_weight"),
        ("gas", "grackle_pressure"),
        ("gas", "grackle_temperature"),
        ("gas", "grackle_dust_temperature"),
    ]
    for field in fields:
        print (f"{field}: {sp[field]}")

    data = {field: sp[field] for field in fields}
    yt.save_as_dataset(ds, filename=f"{output_name}.h5",
                       data=data, extra_attrs=extra_attrs)
