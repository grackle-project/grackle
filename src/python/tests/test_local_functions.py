########################################################################
#
# Grackle python fluid container
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
import json
import numpy as np
import os

from numpy.testing import assert_approx_equal

from gracklepy import \
    FluidContainer, \
    chemistry_data, \
    setup_fluid_container

from gracklepy.utilities.physical_constants import \
    cm_per_mpc, \
    mass_hydrogen_cgs, \
    sec_per_Myr

from testing_common import grackle_data_dir

local_function_test_format_version = 1
_meta_data = {"format_version": local_function_test_format_version}

parameter_grid = {
    "use_grackle": (1,),
    "primordial_chemistry": (0, 1, 2, 3),
    "dust_chemistry": (0, 1),
    "metal_cooling": (0, 1),
    "UVbackground": (0, 1),
    "cmb_temperature_floor": (0, 1),
    "grackle_data_file": ("CloudyData_UVB=HM2012.h5",),
}

exclude_sets = (
    {"dust_chemistry": 1, "metal_cooling": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 0},
    {"dust_chemistry": 1, "primordial_chemistry": 1},
)

grackle_functions = (
    "cooling_time",
    "dust_temperature",
    "gamma",
    "pressure",
    "temperature"
)

function_requirements = {
    "dust_temperature": {"dust_chemistry": 1},
}

base_units = {
    "comoving_coordinates": 1,
    "a_units": 0.1,
    "density_units": mass_hydrogen_cgs,
    "length_units": cm_per_mpc,
    "time_units": sec_per_Myr
}

random_inputs = {
    "density":     {"min": -27, "max":  -20, "log": True},
    "temperature": {"min":   1, "max":    8, "log": True},
    "metallicity": {"min":  -6, "max": -1.5, "log": True},
    "redshift":    {"min":   0, "max":    6, "log": False},
}

def failure_str(my_pars, my_units, my_in, my_out, comp_out, field):
    msg = "Failed local function test:\n"
    msg += "Non-default parameters:\n"
    for par, val in my_pars.items():
        msg += f"\t{par} = {val}\n"
    msg += "Units:\n"
    for unit, val in my_units.items():
        msg += f"\t{unit} = {val}\n"
    msg += "Inputs:\n"
    for my_input, val in my_in.items():
        msg += f"\t{my_input} = {val}\n"
    msg += f"field: {field}\n"
    msg += f"generated output: {my_out[field]}\n"
    msg += f"stored output: {comp_out[field]}"
    return msg

class TestLocalFunctions:
    """
    Tests for the local functions.
    """

    n_input_sets = 10
    seed = 21062024
    digits = 12
    test_file_basename = "local_function_tests.json"

    def setup_parameter_sets(self, generate_answers, test_file):
        """
        Setup parameter-sets for local function tests.

        If we are generating results, then create a list of parameter values
        from the parameter_grid and exclude_sets structures defined above.
        Also, use the base_units dict to create units.

        If we are not generating results, then read everything in from the
        json test file.
        """
        # if we need to reuse this functionality, we need should remove this
        # from the class and convert it into a pytest fixture

        if generate_answers:
            my_sets = self.generate_parameter_sets()
            self.test_sets = [{"parameters": my_set, "units": base_units}
                              for my_set in my_sets]
        else:
            if not os.path.exists(test_file):
                self.skipTest(f"Test file not found: {test_file}.")
            with open(test_file, mode="r") as f:
                self.test_sets = json.load(f)
            load_meta = self.test_sets.pop(0)
            assert load_meta["format_version"] == _meta_data["format_version"], \
              f"Test version mismatch: data file is {load_meta['format_version']}, " + \
              f"source code is {_meta_data['format_version']}."

    def finish_tests(self, generate_answers, test_file):
        """
        Write json test file if we are generating resuls.
        """
        if generate_answers:
            # add the format version to the output
            self.test_sets.insert(0, _meta_data)
            with open(test_file, mode="w") as f:
                f.write(json.dumps(self.test_sets, indent=4))

    def generate_parameter_sets(self):
        """
        Generate all permutations of parameter values except exclusions.
        """

        par_sets = []

        pars = parameter_grid.keys()
        for perm in itertools.product(*parameter_grid.values()):
            my_dict = dict(zip(pars, perm))

            exclude = False
            for my_exclude in exclude_sets:
                matches = sum([my_exclude[par] == my_dict[par]
                               for par in my_exclude])
                if matches == len(my_exclude):
                    exclude = True
                    break

            if exclude:
                continue

            par_sets.append(my_dict)

        return par_sets

    def generate_base_inputs(self):
        """
        Generate sets of sensible random inputs using random_inputs dict.

        This generates just the minimal inputs, out of which a
        full set of input fields can be generated with a given
        set of runtime parameters.
        """

        values = {}
        rng = np.random.default_rng(self.seed)
        for field, config in random_inputs.items():
            vals = (config["max"] - config["min"]) * \
              rng.random(size=self.n_input_sets) + config["min"]
            if config["log"]:
                vals = 10**vals
            values[field] = vals

        base_inputs = [{field: values[field][i] for field in values}
                       for i in range(self.n_input_sets)]
        return base_inputs

    def test_local_functions(self, answertestspec):
        test_file = os.path.join(answertestspec.answer_dir, self.test_file_basename)
        self.setup_parameter_sets(answertestspec.generate_answers, test_file)
        for test_set in self.test_sets:

            par_set = test_set["parameters"]
            base_units = test_set["units"]

            my_chemistry = chemistry_data()
            for par, val in par_set.items():
                setattr(my_chemistry, par, val)

            # add path to data file here so we only have to save the base name
            if "grackle_data_file" in par_set:
                setattr(my_chemistry, "grackle_data_file",
                        os.path.join(grackle_data_dir, par_set["grackle_data_file"]))

            if answertestspec.generate_answers:
                my_tests = []
                base_inputs = self.generate_base_inputs()
                n_input_sets = len(base_inputs)
            else:
                my_tests = test_set["tests"]
                n_input_sets = len(my_tests)
            assert n_input_sets == self.n_input_sets

            for itest in range(n_input_sets):
                my_units = base_units.copy()
                if answertestspec.generate_answers:
                    base_input_set = base_inputs[itest]
                    redshift = base_input_set["redshift"]
                else:
                    my_in = my_tests[itest]["input"]
                    redshift = my_in.pop("redshift")

                my_units["a_value"] = 1 / (1 + redshift) / my_units["a_units"]
                for unit, val in my_units.items():
                    setattr(my_chemistry, unit, val)
                my_chemistry.set_velocity_units()

                if answertestspec.generate_answers:
                    fc = setup_fluid_container(
                        my_chemistry,
                        density=base_input_set["density"],
                        temperature=base_input_set["temperature"],
                        metal_mass_fraction=base_input_set["metallicity"],
                        converge=False)
                    my_in = {field: fc[field][0] for field in fc.input_fields}
                    my_in["redshift"] = redshift

                else:
                    rval = my_chemistry.initialize()
                    if rval == 0:
                        raise RuntimeError("Failed to initialize chemistry_data.")
                    fc = FluidContainer(my_chemistry, 1)
                    for field in my_in:
                        fc[field][:] = my_in[field]

                # Call local functions and store results.
                my_out = {}
                for fname in grackle_functions:
                    my_reqs = function_requirements.get(fname, {})
                    reqs_met = sum([my_reqs[p] == getattr(my_chemistry, p, None)
                                   for p in my_reqs])
                    if reqs_met != len(my_reqs):
                        continue
                    getattr(fc, f"calculate_{fname}")()
                    my_out[fname] = fc[fname][0]

                # Compare with existing results unless we are generating them.
                if answertestspec.generate_answers:
                    my_tests.append({"input": my_in, "output": my_out})
                else:
                    comp_out = my_tests[itest]["output"]
                    for field in my_out:
                        err_msg=failure_str(
                            par_set, my_units, my_in,
                            my_out, comp_out, field)
                        assert_approx_equal(
                            comp_out[field], my_out[field],
                            significant=self.digits,
                            err_msg=err_msg)

            if answertestspec.generate_answers:
                test_set["tests"] = my_tests
        self.finish_tests(answertestspec.generate_answers, test_file)
