
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

import os
import pytest
import subprocess

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

code_examples = ["c_example",
                 "c_local_example",
                 "cxx_example",
                 "cxx_omp_example",
                 "cxx_grid_example",
                 "fortran_example"]


def run_command(command, cwd, env):
    try:
        subprocess.check_output(
            command.split(' '), stderr=subprocess.STDOUT,
            cwd=cwd, env=env)
    except subprocess.CalledProcessError as er:
        raise RuntimeError(
            f"Command {command} failed with return code {er.returncode} "
            f"and the following output: {er.output}")

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

    dflt_choice = "classic" if _USING_TRADITIONAL_BUILD == True else ""
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
    run_command(command, examples_path, env)

    # test that example compiles
    assert os.path.exists(os.path.join(examples_path, example))

    # try to run the example code
    command = f"./{example}"
    run_command(command, examples_path, env)

    command = f"{make_command} clean"
    run_command(command, examples_path, env)
