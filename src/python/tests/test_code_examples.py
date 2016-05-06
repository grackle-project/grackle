
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
import subprocess

current_path = os.path.abspath(__file__)

examples_path = os.path.join(os.path.dirname(current_path), "../..", "example")

code_examples = ["c_example",
                 "cxx_example",
                 "cxx_omp_example"]


def run_command(command, cwd, env):
    try:
        subprocess.check_output(
            command.split(' '), stderr=subprocess.STDOUT,
            cwd=cwd, env=env)
    except subprocess.CalledProcessError as er:
        raise RuntimeError('Command %s failed with return code %s '
                           'and the following output: %s' %
                           (command, er.returncode, er.output))

def code_example_test(example):
    env = dict(os.environ)
    curdir = os.getcwd()
    os.chdir(examples_path)
    command = 'make %s' % example
    run_command(command, examples_path, env)

    # test that example compiles
    assert os.path.exists(example)

    # try to run the example code
    command = "./%s" % example
    run_command(command, examples_path, env)

    command = "make clean"
    run_command(command, examples_path, env)
    
    os.chdir(curdir)

def test_code_examples():
    for example in code_examples:
        yield code_example_test, example
