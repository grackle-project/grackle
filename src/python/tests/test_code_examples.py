
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

    _build_dir = os.getenv("PYTEST_CMAKE_BUILD_DIR", "")
    if (_USING_TRADITIONAL_BUILD in [None, True]) and _build_dir == "":
        # lets assume that the underlying library was constructed with the
        # traditional build-system
        make_command = "make"
    elif _USING_TRADITIONAL_BUILD and _build_dir != "":
        raise RuntimeError("PYTEST_CMAKE_BUILD_DIR env variable is set when "
                           "libgrackle is built with traditional build system")
    elif (_USING_TRADITIONAL_BUILD is False) and _build_dir == "":
        raise RuntimeError("PYTEST_CMAKE_BUILD_DIR env variable is needed when "
                           "libgrackle is built with cmake")
    else:
        # this branch is executed if the underlying libgrackle library was
        # instead constructed with CMake
        build_dir = os.path.join(os.getcwd(), _build_dir)
        make_command = (
            f"make -f Makefile.out-of-source CMAKE_BUILD_DIR={build_dir}")

    env = dict(os.environ)
    curdir = os.getcwd()
    os.chdir(examples_path)
    command = f'{make_command} {example}'
    run_command(command, examples_path, env)

    # test that example compiles
    assert os.path.exists(example)

    # try to run the example code
    command = f"./{example}"
    run_command(command, examples_path, env)

    command = f"{make_command} clean"
    run_command(command, examples_path, env)
    
    os.chdir(curdir)
