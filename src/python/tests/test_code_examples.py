
########################################################################
#
# Tests the code examples
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import json
import os
import pytest
import re
import subprocess

from pygrackle.utilities.testing import \
    generate_test_results, \
    test_answers_dir

try:
    from pygrackle.grackle_wrapper import uses_in_source_build
    _USING_TRADITIONAL_BUILD = uses_in_source_build()
except ImportError:
    # this branch is executed when grackle was compiled with OpenMP
    # -> last time I checked, this won't happen if OpenMP is used in a cmake,
    #    build. But we won't depend on that behavior
    _USING_TRADITIONAL_BUILD = None

current_path = os.path.abspath(__file__)
current_dir = os.path.dirname(current_path)

examples_path = os.path.join(current_dir, "../..", "example")

code_examples = (
    "c_example",
    "c_local_example",
    "cxx_example",
    "cxx_omp_example",
    "cxx_grid_example",
    "fortran_example",
)

compare_exclude = (
    "cxx_omp_example",
    "cxx_grid_example"
)

rfields = (
    "cooling_time",
    "dust_temperature",
    "gamma",
    "pressure",
    "temperature"
)

test_file = os.path.join(test_answers_dir, "code_examples.json")
if generate_test_results and os.path.exists(test_file):
    os.remove(test_file)
if not generate_test_results and not os.path.exists(test_file):
    raise RuntimeError(
        f"Code example results file not found: {test_file}")

def run_command(command, cwd, env, timeout=None):
    proc = subprocess.run(
        command, shell=True,
        cwd=cwd, env=env, capture_output=True)
    if proc.returncode == 0:
        return proc
    else:
        raise RuntimeError(
            f"Command {command} failed with return value "
            f"{proc.returncode} and the following stderr output "
            f"{proc.stderr}")

def parse_output(ostr):
    results = {field: None for field in rfields}

    if isinstance(ostr, bytes):
        ostr = ostr.decode("utf8")
    lines = ostr.split("\n")
    for line in lines:
        for field in results:
            match = re.match(f"^ ?{field} = ", line)
            if match is None:
                continue
            _, rside = line.split(" = ")
            rparts = rside.split()
            if len(rparts) == 1:
                val = rparts[0][:-1]
            elif len(rparts) == 2:
                val = rparts[0]
            else:
                raise RuntimeError(
                    f"Cannot grab field values from line: {line}.")

            if results[field] is not None:
                raise RuntimeError(
                    f"Already have value for {field}.")
            results[field] = val

    return results

@pytest.mark.parametrize("example", code_examples)
def test_code_examples(example):
    # under the classic build system, we could just execute `make` in the examples
    # directory and there is a good chance that it would work out...
    # -> now we require the PYTEST_CODE_LINK_CHOICE environment variable to be
    #    explicitly set in order to dictate how to execute the test
    # -> if the variable isn't set, then we skip the test
    # -> purely for backwards compatability, when the PYTEST_CODE_LINK_CHOICE
    #    variable isn't specified, but pygrackle was built with the classic system
    #    then we assume that the user wants to test stuff built with the classic
    #    build-system... (Hopefully, we can remove this extra logic in future
    #    revisions -- it would simplify a lot!)

    dflt_choice = "classic" if _USING_TRADITIONAL_BUILD else ""
    choice = os.getenv("PYTEST_CODE_LINK_CHOICE", dflt_choice)

    if choice == "":
        pytest.skip("the 'PYTEST_CODE_LINK_CHOICE' environment variable must be "
                    "set to run this test")
    elif choice == "classic":
        make_command = "make"
    elif choice.startswith("cmake:") and (len(choice) > 6):
        build_dir = os.path.expanduser(choice[6:])
        if not os.path.isabs(build_dir):
            build_dir = os.path.join(os.getcwd(), build_dir)
        if not os.path.isdir(build_dir):
            raise RuntimeError(f"{build_dir} specified as path to cmake-build dir, "
                               "but it doesn't exist")
        make_command = (
            f"make -f Makefile.out-of-source CMAKE_BUILD_DIR={build_dir}")
    else:
        raise RuntimeError("PYTEST_CODE_LINK_CHOICE must be '', 'classic' or "
                           "'cmake:<path/to/build>'. {choice!r} is invalid")
    env = dict(os.environ)
    command = f'{make_command} {example}'
    run_command(command, examples_path, env, timeout=60)

    # test that example compiles
    assert os.path.exists(os.path.join(examples_path, example))

    # try to run the example code
    command = f"./{example}"
    proc = run_command(command, examples_path, env, timeout=120)
    if example not in compare_exclude:
        results = parse_output(proc.stdout)

        if generate_test_results:
            if os.path.exists(test_file):
                with open(test_file, mode="r") as f:
                    all_results = json.load(f)
            else:
                all_results = {}

            all_results.update({example: results})
            with open(test_file, mode="w") as f:
                f.write(json.dumps(all_results, indent=4))

        else:
            with open(test_file, mode="r") as f:
                all_results = json.load(f)

            comp_results = all_results[example]
            for field in comp_results:
                err_msg = f"In {example}: mismatch for {field} - " + \
                  f"old: {comp_results[field]}, new: {results[field]}"
                assert comp_results[field] == results[field], err_msg

    command = f"{make_command} clean"
    run_command(command, examples_path, env, timeout=60)
