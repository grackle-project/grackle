
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

from testing_common import grackle_install_dir

try:
    from pygrackle.__config__ import \
        _GrackleBuild, \
        _grackle_build
    _USING_TRADITIONAL_BUILD = _grackle_build == _GrackleBuild.ExternalClassic
except ImportError:
    # this branch is executed when grackle was compiled with OpenMP
    # -> last time I checked, this won't happen if OpenMP is used in a cmake,
    #    build. But we won't depend on that behavior
    _USING_TRADITIONAL_BUILD = None

examples_dir = os.path.join(grackle_install_dir, "src", "example")

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

@pytest.fixture(scope="module")
def test_file(answertestspec):
    # this fixture has module-scope so that we will delete the test file once
    test_file = os.path.join(answertestspec.answer_dir, "code_examples.json")
    if answertestspec.generate_answers and os.path.exists(test_file):
        os.remove(test_file)

    # if test_file doesn't exist and answertestspec.generate_answers is False,
    # defer any error reporting until within the test-case (after the test-case
    # determines whether or not it should be skipped)
    return test_file

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
    results = {when: {field: None for field in rfields}
               for when in ("Before", "After")}

    if isinstance(ostr, bytes):
        ostr = ostr.decode("utf8")
    lines = ostr.split("\n")
    for line in lines:
        match = re.match(r"^ ?(\w+) - (\w+) = ", line)
        if match is None:
            continue
        when, field = match.groups()
        _, rside = line.split(" = ")
        rparts = rside.split()
        if len(rparts) == 1:
            val = rparts[0][:-1]
        elif len(rparts) == 2:
            val = rparts[0]
        else:
            raise RuntimeError(
                f"Cannot grab field values from line: {line}.")

        if results[when][field] is not None:
            raise RuntimeError(
                f"Already have value for {field}.")
        results[when][field] = val

    return results

@pytest.mark.parametrize("example", code_examples)
def test_code_examples(answertestspec, test_file, example):
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

    # if we aren't generating test-answers, and the test-file can't be found
    # report an error (we explicitly wait to do this until after we have
    # decided whether to skip the test or not).
    if not answertestspec.generate_answers and not os.path.exists(test_file):
        raise RuntimeError(f"Code example results file not found: {test_file}")

    env = dict(os.environ)

    # compile the example
    command = f'{make_command} {example}'
    run_command(command, examples_dir, env, timeout=60)

    # test that example compiles
    assert os.path.exists(os.path.join(examples_dir, example))

    # try to run the example code
    command = f"./{example}"
    proc = run_command(command, examples_dir, env, timeout=120)
    if example not in compare_exclude:
        results = parse_output(proc.stdout)

        if answertestspec.generate_answers:
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
            failures = 0
            err_msg = f"{example}:\n"
            for when in comp_results:
                for field in comp_results[when]:
                    if comp_results[when][field] == results[when][field]:
                        continue
                    failures += 1
                    err_msg += f"\t{when} - {field} - " + \
                      f"old: {comp_results[when][field]}, " + \
                      f"new: {results[when][field]}\n"
            assert failures == 0, err_msg

    command = f"{make_command} clean"
    run_command(command, examples_dir, env, timeout=60)
