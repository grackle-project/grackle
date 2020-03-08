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

import numpy as np
import os
import yt

from pygrackle import \
    add_grackle_fields

DS_NAME = "IsolatedGalaxy/galaxy0030/galaxy0030"

if 'YT_DATA_DIR' in os.environ:
    ds_path = os.sep.join([os.environ['YT_DATA_DIR'], DS_NAME])
else:
    ds_path = DS_NAME

ds = yt.load(ds_path)

my_dir = os.path.dirname(os.path.abspath(__file__))
grackle_data_file = bytes(os.path.join(
    my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012.h5"), 'utf-8')

grackle_pars = {'grackle_data_file': grackle_data_file,
                'UVbackground': 1}

add_grackle_fields(ds, parameters=grackle_pars)

sp = ds.sphere(ds.domain_center, (10, 'kpc'))
print (sp['gas', 'grackle_cooling_time'])
print (sp['gas', 'grackle_gamma'])
print (sp['gas', 'grackle_pressure'])
print (sp['gas', 'grackle_temperature'])
print (sp['gas', 'grackle_dust_temperature'])
