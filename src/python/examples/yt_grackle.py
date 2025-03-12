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
from pygrackle.utilities.model_tests import parse_model_cliargs

output_name = os.path.basename(__file__[:-3]) # strip off ".py"

DS_NAME = "IsolatedGalaxy/galaxy0030/galaxy0030"

if 'YT_DATA_DIR' in os.environ:
    ds_path = os.sep.join([os.environ['YT_DATA_DIR'], DS_NAME])
else:
    ds_path = DS_NAME

if __name__ == "__main__":
    # If we are running the script through the testing framework, then we will
    # pass in arguments (recall sys.argv is just the file name)
    if len(sys.argv) > 1:
        tmp = parse_model_cliargs(sys.argv[1:], output_name)
        if tmp is None:
            sys.exit(0) # exit early
        output_name = tmp.output_name
        extra_attrs = tmp.extra_attrs

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
