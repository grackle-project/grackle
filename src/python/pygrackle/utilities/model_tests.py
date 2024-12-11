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
import pprint
from typing import NamedTuple

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

def _build_base_parser():
    parser = argparse.ArgumentParser(
        description=(
            "This is a script that provides an illustrative example of how to "
            "use Grackle can be used in a simplified model. When using the "
            "script for learning, you should invoke it without any arguments. "
            "Some subcommands are implemented to support using this script "
            "for testing purposes"
        )
    )

    subparsers = parser.add_subparsers(
        required=True,
        description = (
            "these sub-commands support the use of the script's logic as part "
            "of the test suite. We implement a command for evaluating the "
            "logic with a parameter presets and an input-grid preset (both "
            "sets of presets are associated with integers). Subcommands are "
            "also provided to query the presets."
        )
    )

    parser_run_test = subparsers.add_parser(
        "run-test", description="execute the logic as part of a test-case"
    )
    parser_run_test.set_defaults(kind="run-test")
    _kw = {'type' : int, 'required' : True}
    parser_run_test.add_argument(
        "--param-preset", help="index of the parameter preset", **_kw
    )
    parser_run_test.add_argument(
        "--input-preset", help="index of the input-grid preset", **_kw
    )
    parser_run_test.add_argument(
        "--skip-plot", action="store_true", help="skip any plotting"
    )

    _triples = [
        ("show-param-preset", "show the parameter preset", True),
        ("show-input-preset", "show the input-grid preset", True),
        ("num-param-preset", "number of parameter presets", False),
        ("num-input-preset", "number of input-grid presets", False),
    ]
    for name, descr, has_arg in _triples:
        p = subparsers.add_parser(name, description=descr)
        p.set_defaults(kind=name)
        if has_arg:
            p.add_argument("index", type=int, help="index of the preset")
    return parser


_MODEL_CLI_PARSER = _build_base_parser()

class ParsedModelProps(NamedTuple):
    my_chemistry: chemistry_data
    input_set: dict
    skip_plot: bool
    output_name: str
    extra_attrs: dict

def _get_preset(index_value, sets, index_name):
    try:
        return sets[index_value]
    except IndexError:
        if len(sets) == 1:
            raise ValueError(f"{index_name} must be 0") from None
        raise ValueError(
            f"{index_name} must satisfy 0 <= index < {len(sets)}"
        ) from None

def parse_model_cliargs(arg_l, model_name):
    """
    Parses command line arguments and takes appropriate actions

    Parameters
    ----------
    arg_l: list of str
        The list of cli args
    model_name: str
        The name of the model to be tested

    Returns
    -------
    out : None or ParsedModelProps
        If a subcommand was invoked that implies that the program should
        immediately exit, this returns None. Otherwise, this provides the
        appropriate data to be used while invoking the simple model's logic

    Notes
    -----
    Rather than defining the _model_test_grids variable inside of this file and
    matching up the keys with specified model_name argument, I think it might
    be more robust for the possible _model_test_grids to be directly passed
    into this function.
    """
    if model_name not in model_sets:
        raise ValueError(f"Unkown model name: {model_name}.")
    my_model = model_sets[model_name]
    par_sets, input_sets = my_model["parameter_sets"], my_model["input_sets"]

    args = _MODEL_CLI_PARSER.parse_args(arg_l)

    if args.kind == "num-param-preset":
        print(len(par_sets))
        return None
    elif args.kind == "num-input-preset":
        print(len(input_sets))
        return None
    elif args.kind == "show-param-preset":
        par_set = _get_preset(args.index, par_sets, 'preset index')
        pprint.pprint(par_set)
        return None
    elif args.kind == "show-input-preset":
        input_set = _get_preset(args.index, input_sets, 'preset index')
        pprint.pprint(input_set)
        return None
    elif args.kind == "run-test":
        par_index = args.param_preset
        par_set = _get_preset(par_index, par_sets, '--param-preset')
        my_chemistry = chemistry_data()
        for par, val in par_set.items():
            if par == "grackle_data_file":
                val = os.path.join(grackle_data_dir, val)
            setattr(my_chemistry, par, val)

        input_index = args.input_preset
        input_set = _get_preset(input_index, input_sets, '--input-preset')
        return ParsedModelProps(
            my_chemistry=my_chemistry,
            input_set=input_set,
            skip_plot=args.skip_plot,
            output_name=f"{model_name}_{par_index}_{input_index}",
            extra_attrs={"format_version": model_test_format_version}
        )
    else:
        raise RuntimeError("something unexepected happened!")
