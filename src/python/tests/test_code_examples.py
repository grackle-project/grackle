
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

current_path = os.path.abspath(__file__)

examples_path = os.path.join(os.path.dirname(current_path), "../..", "example")

code_examples = ["c_example",
                 "c_local_example",
                 "cxx_example",
                 "cxx_omp_example",
                 "cxx_grid_example",
                 "fortran_example",
                 "test_calc_temp1d_cloudy",
                 "test_interpolators"]


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
    env = dict(os.environ)
    curdir = os.getcwd()
    os.chdir(examples_path)
    command = f'make {example}'
    run_command(command, examples_path, env)

    # test that example compiles
    assert os.path.exists(example)

    # try to run the example code
    command = f"./{example}"
    run_command(command, examples_path, env)

    command = "make clean"
    run_command(command, examples_path, env)
    
    os.chdir(curdir)
