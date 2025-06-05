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

import argparse
import itertools
import os
import sys

from pygrackle import chemistry_data

if sys.version_info < (3, 9):
    from functools import lru_cache as cache
else:
    from functools import cache

model_test_format_version = 1

# These are generally applicable to all scripts.
_parameter_exclude = (
    {"dust_chemistry": 1, "metal_cooling": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 1},
    {"metal_chemistry": 1, "primordial_chemistry": 0},
    {"metal_chemistry": 1, "primordial_chemistry": 1},
    {"metal_chemistry": 1, "metal_cooling": 0},
    {"grackle_data_file": "CloudyData_noUVB.h5", "UVbackground": 1},
)

_model_test_grids = \
{
    "cooling_rate": \
    {
        "standard_variants": \
        {
            "parameters":
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "primordial_chemistry": 4,
                    "dust_chemistry": 1,
                    "metal_cooling": 1,
                    "UVbackground": 1,
                    "cmb_temperature_floor": 1,
                    "grackle_data_file": "CloudyData_UVB=HM2012.h5",
                },
                "variants": \
                {
                    "UVbackground": (0,),
                    "cmb_temperature_floor": (0,),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 1.,
                    "redshift": 0.,
                    "specific_heating_rate": 0.,
                    "volumetric_heating_rate": 0.,
                },
                "variants": \
                {
                    "metallicity": (0.,),
                    "redshift": (2.,),
                    "specific_heating_rate": (1.,),
                    "volumetric_heating_rate": (1e-24,),
                },
            }
        },

        "standard_combinations": \
        {
            "parameters": \
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "metal_cooling": 1,
                    "UVbackground": 1,
                    "cmb_temperature_floor": 1,
                    "grackle_data_file": "CloudyData_UVB=HM2012.h5",
                },
                "grid": \
                {
                    "primordial_chemistry": (0, 1, 2, 3, 4),
                    "dust_chemistry": (0, 1),
                    "metal_cooling": (0, 1),
                    "metal_chemistry": (0, 1),
                },
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 1.,
                    "redshift": 0.,
                    "specific_heating_rate": 0.,
                    "volumetric_heating_rate": 0.,
                }
            }
        },

        "cloudy_table_combinations": \
        {
            "parameters":
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "metal_cooling": 1,
                    "cmb_temperature_floor": 1,
                },
                "grid": \
                {
                    "primordial_chemistry": (0, 4),
                    "UVbackground": (0, 1),
                    "grackle_data_file": ("CloudyData_noUVB.h5",
                                          "CloudyData_UVB=FG2011_shielded.h5",
                                          "CloudyData_UVB=FG2011.h5",
                                          "CloudyData_UVB=HM2012_high_density.h5",
                                          "CloudyData_UVB=HM2012_shielded.h5",
                                          "CloudyData_UVB=HM2012.h5"),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 1.,
                    "redshift": 0.,
                    "specific_heating_rate": 0.,
                    "volumetric_heating_rate": 0.,
                },
            }
        },

        "self_shielding_combinations": \
        {
            "parameters":
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "metal_cooling": 1,
                    "cmb_temperature_floor": 1,
                    "UVbackground": 1,
                },
                "grid": \
                {
                    "primordial_chemistry": (0, 4),
                    "grackle_data_file": ("CloudyData_UVB=FG2011_shielded.h5",
                                          "CloudyData_UVB=HM2012_shielded.h5"),
                    "self_shielding_method": (1, 2, 3),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 0.1,
                    "redshift": 2.,
                    "specific_heating_rate": 0.,
                    "volumetric_heating_rate": 0.,
                },
            }
        },
    },

    "cooling_cell": \
    {
        "standard_combinations": \
        {
            "parameters": \
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "metal_cooling": 1,
                    "cmb_temperature_floor": 1,
                    "grackle_data_file": "CloudyData_UVB=HM2012.h5",
                    "temperature_floor_scalar": 5e4,
                },
                "grid": \
                {
                    "primordial_chemistry": (0, 1, 2, 3),
                    "metal_chemistry": (0, 1),
                    "UVbackground": (0, 1),
                    "use_temperature_floor": (0, 1),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 0.1,
                    "redshift": 0.,
                }
            }
        }
    },

    "freefall": \
    {
        "standard_combinations": \
        {
            "parameters": \
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "photoelectric_heating": 0,
                    "dust_recombination_cooling": 0,
                    "cmb_temperature_floor": 1,
                    "CaseBRecombination": 1,
                    "cie_cooling": 1,
                    "h2_optical_depth_approximation": 1,
                    "use_primordial_continuum_opacity": 1,
                    "grackle_data_file": "cloudy_metals_2008_3D.h5",
                },
                "grid": \
                {
                    "primordial_chemistry": (2, 3, 4),
                    "metal_cooling": (0, 1),
                    "dust_chemistry": (0, 1),
                }
            },
            "inputs": \
            {
                "grid": \
                {
                    "metallicity": (1e-10, 1e-4, 1e-3),
                }
            }
        },

        "primordial_rate_variants": \
        {
            "parameters": \
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "primordial_chemistry": 4,
                    "CaseBRecombination": 1,
                    "cie_cooling": 1,
                    "h2_optical_depth_approximation": 1,
                    "use_primordial_continuum_opacity": 1,
                    "h2_charge_exchange_rate": 1,
                    "h2_h_cooling_rate": 1,
                    "h2_cooling_rate": 2,
                    "hd_reaction_rates": 0,
                    "H2_self_shielding": 0,
                    "three_body_rate": 0,
                },
                "variants": \
                {
                    "CaseBRecombination": (0,),
                    "cie_cooling": (0, 2),
                    "h2_optical_depth_approximation": (0,),
                    "use_primordial_continuum_opacity": (0,),
                    "h2_charge_exchange_rate": (2,),
                    "h2_cooling_rate": (0, 1, 3),
                    "h2_h_cooling_rate": (2,),
                    "hd_reaction_rates": (1,),
                    "H2_self_shielding": (3,),
                    "three_body_rate": (1, 2, 3, 4, 5),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 0.,
                }
            }
        },

        "metal_dust_chemistry_variants": \
        {
            "parameters": \
            {
                "defaults": \
                {
                    "use_grackle": 1,
                    "primordial_chemistry": 4,
                    "CaseBRecombination": 1,
                    "cie_cooling": 2,
                    "h2_optical_depth_approximation": 1,
                    "use_primordial_continuum_opacity": 1,
                    "h2_charge_exchange_rate": 1,
                    "h2_h_cooling_rate": 1,
                    "h2_cooling_rate": 3,
                    "hd_reaction_rates": 1,
                    "grackle_data_file": "cloudy_metals_2008_3D.h5",
                    "metal_cooling": 1,
                    "metal_chemistry": 1,
                    "three_body_rate": 4,
                    "dust_chemistry": 1,
                    "use_dust_density_field": 1,
                    "photoelectric_heating": 0,
                    "dust_recombination_cooling": 0,
                    "dust_sublimation": 1,
                    "grain_growth": 1,
                    "dust_species": 1,
                    "use_multiple_dust_temperatures": 1,
                    "multi_metals": 0,
                    "metal_abundances": 0,
                },
                "variants": \
                {
                    "dust_species": (2, 3),
                    "multi_metals": (1,),
                    "metal_abundances": (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11),
                }
            },
            "inputs": \
            {
                "defaults": \
                {
                    "metallicity": 1e-3,
                }
            }
        }

    },

    "yt_grackle": \
    {
        "standard": \
        {
            "parameters": {},
            "inputs": {}
        }
    }
}

def generate_value_sets(config, exclude_sets=None):
    """
    Generate all sets of values excepting exclusions.

    These sets can be generated either as all permutations of a grid
    or as a default and a set of single-parameter variations.
    """

    if exclude_sets is None:
        exclude_sets = {}

    par_defaults = config.get("defaults", {})

    if "grid" in config and "variants" in config:
        raise ValueError("Must provide grid or variants, but not both.")

    my_sets = []

    # do all permutations of a grid of parameters
    if "grid" in config:
        pars = config["grid"].keys()
        for perm in itertools.product(*config["grid"].values()):
            my_sets.append(dict(zip(pars, perm)))

    # do variants of the default parameters one at a time
    elif "variants" in config:
        # do the default
        my_sets.append({})

        # add the variants
        for par, vals in config["variants"].items():
            for val in vals:
                my_sets.append({par: val})

    # just do the defaults
    else:
        my_sets.append({})

    # remove excluded parameter combinations
    par_sets = []
    for my_set in my_sets:
        my_dict = par_defaults.copy()
        my_dict.update(my_set)

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

def generate_model_sets():
    model_store = {}
    model_parametrization = []

    for model_name, model_variants in _model_test_grids.items():
        model_store[model_name] = {}
        for model_label, model in model_variants.items():
            my_model = {}
            if "parameter_sets" not in model:
                my_model["parameter_sets"] = \
                  generate_value_sets(model["parameters"],
                                      exclude_sets=_parameter_exclude)
            if "input_sets" not in model:
                my_model["input_sets"] = \
                  generate_value_sets(model["inputs"])
            model_store[model_name][model_label] = my_model

            for par_index in range(len(my_model["parameter_sets"])):
                for input_index in range(len(my_model["input_sets"])):
                    model_parametrization.append(
                        (model_name, model_label, par_index, input_index))

    return model_store, model_parametrization


model_sets, model_parametrization = generate_model_sets()

def get_model_set(model_name, model_variant, parameter_index, input_index):
    """
    Create objects and variables to be used in testing one
    of the Python example scripts.
    """

    if model_name not in model_sets:
        raise ValueError(f"Unknown model name: {model_name}.")

    if model_variant not in model_sets[model_name]:
        raise ValueError(
            f"Unknown model variant in {model_name}: {model_variant}.")

    my_model = model_sets[model_name][model_variant]
    par_set = my_model["parameter_sets"][parameter_index]
    input_set = my_model["input_sets"][input_index]

    return par_set, input_set

@cache
def _build_parser():
    p = argparse.ArgumentParser(
        description=(
            "This script illustrates an example of using Grackle. To use the "
            "the script for learning, invoke it without any args (the command "
            "line arguments are ONLY for tests)"
        )
    )
    p.add_argument("--out-dir", default=".", help="the output directory")
    p.add_argument("model_variant", type=str)
    p.add_argument("parameter_index", type=int)
    p.add_argument("input_index", type=int)
    return p

def get_test_variables(model_name, arg_l):
    """
    Setup objects and variables for a model test.

    This is called from one of the Python example scripts
    when being run from pytest. This avoids cluttering the
    script with a bunch of code that a user just wanting to
    play with it shouldn't have to see.
    """

    # we import this here, rather than at global scope to let us import the
    # module when we don't have an editable install (for testing purposes)
    from pygrackle.utilities.data_path import grackle_data_dir

    args = _build_parser().parse_args(args=arg_l)
    model_variant = args.model_variant
    par_index, input_index = args.parameter_index, args.input_index

    par_set, input_set = get_model_set(
        model_name, model_variant, par_index, input_index)

    # setup chemistry data object
    my_chemistry = chemistry_data()
    for par, val in par_set.items():
        if par == "grackle_data_file":
            val = os.path.join(grackle_data_dir, val)
        setattr(my_chemistry, par, val)

    output_name = os.path.join(
        args.out_dir, f"{model_name}_{model_variant}_{par_index}_{input_index}"
    )
    extra_attrs = {"format_version": model_test_format_version}
    extra_attrs.update(par_set)
    extra_attrs.update(input_set)
    extra_attrs = {k: str(v) for k, v in extra_attrs.items()}

    # all variables to be set as globals in the script
    my_globals = {
        "my_chemistry": my_chemistry,
        "output_name": output_name,
        "extra_attrs": extra_attrs,
    }
    # input variables used in the script
    my_globals.update(input_set)

    return my_globals
