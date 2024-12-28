########################################################################
#
# Infrastructure for testing Python example scripts
#
#
# Copyright (c) Enzo/Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import itertools
import os

from pygrackle import chemistry_data
from pygrackle.utilities.data_path import grackle_data_dir

model_test_format_version = 1

# These are generally applicable to all scripts.
_parameter_exclude = (
    {"dust_chemistry": 1, "metal_cooling": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 1},
)

_model_test_grids = \
{
    "cooling_rate": \
    {
        "parameter_grid": \
        {
            "use_grackle": (1,),
            "primordial_chemistry": (0, 1, 2, 3),
            "dust_chemistry": (1,),
            "metal_cooling": (1,),
            "UVbackground": (0, 1),
            "cmb_temperature_floor": (1,),
            "grackle_data_file": ("CloudyData_UVB=HM2012.h5",),
        },
        "input_grid": \
        {
            "metallicity": (0., 1.),
            "redshift": (0., 2.),
            "specific_heating_rate": (0.,),
            "volumetric_heating_rate": (0.,),
        }
    },
    "cooling_cell": \
    {
        "parameter_grid": \
        {
            "use_grackle": (1,),
            "primordial_chemistry": (0, 1, 2, 3),
            "metal_cooling": (1,),
            "UVbackground": (1,),
            "cmb_temperature_floor": (1,),
            "grackle_data_file": ("CloudyData_UVB=HM2012.h5",),
        },
        "input_grid": \
        {
            "metallicity": (0.1,),
            "redshift": (0.,),
        }
    },
    "freefall": \
    {
        "parameter_grid": \
        {
            "use_grackle": (1,),
            "primordial_chemistry": (2, 3),
            "metal_cooling": (0, 1),
            "dust_chemistry": (0, 1),
            "photoelectric_heating": (0,),
            "cmb_temperature_floor": (1,),
            "CaseBRecombination": (1,),
            "cie_cooling": (1,),
            "h2_optical_depth_approximation": (1,),
            "grackle_data_file": ("cloudy_metals_2008_3D.h5",),
        },
        "input_grid": \
        {
            "metallicity": (1e-10, 1e-4, 1e-3),
        }
    },
    "yt_grackle": \
    {
        "parameter_grid": {},
        "input_grid": {}
    }
}

def generate_value_sets(parameter_grid, exclude_sets=None):
    """
    Generate all permutations of values excepting exclusions.
    """

    if exclude_sets is None:
        exclude_sets = {}

    par_sets = []
    pars = parameter_grid.keys()
    for perm in itertools.product(*parameter_grid.values()):
        my_dict = dict(zip(pars, perm))

        exclude = False
        for my_exclude in exclude_sets:
            matches = sum([my_exclude[par] == my_dict[par]
                           for par in my_exclude if par in my_dict])
            if matches == len(my_exclude):
                exclude = True
                break

        if exclude:
            continue

        par_sets.append(my_dict)
    return par_sets

def generate_model_sets(model_store):
    for model_name, model in _model_test_grids.items():
        my_model = {}
        if "parameter_sets" not in model:
            my_model["parameter_sets"] = \
              generate_value_sets(model["parameter_grid"],
                                  exclude_sets=_parameter_exclude)
        if "input_sets" not in model:
            my_model["input_sets"] = \
              generate_value_sets(model["input_grid"])
        model_store[model_name] = my_model


model_sets = {}
generate_model_sets(model_sets)

def get_model_set(model_name, parameter_index, input_index):
    """
    Create objects and variables to be used in testing one
    of the Python example scripts.
    """

    if model_name not in model_sets:
        raise ValueError("Unkown model name: {model_name}.")

    my_model = model_sets[model_name]
    par_set = my_model["parameter_sets"][parameter_index]
    input_set = my_model["input_sets"][input_index]

    my_chemistry = chemistry_data()
    for par, val in par_set.items():
        if par == "grackle_data_file":
            val = os.path.join(grackle_data_dir, val)
        setattr(my_chemistry, par, val)
    return my_chemistry, input_set
