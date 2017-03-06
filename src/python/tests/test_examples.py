
########################################################################
#
# Tests for python examples
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import contextlib
import glob
import os
import pytest
import shutil
import subprocess
import tempfile
import yt

from pygrackle.utilities.testing import \
    assert_allclose

@contextlib.contextmanager
def temporary_directory():
    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=curdir)
    try:
        yield tmpdir
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)


current_path = os.path.abspath(__file__)

EXAMPLES_GLOB = [os.path.dirname(
    os.path.dirname(current_path)), 'examples', '*.py']

python_examples = glob.glob(os.sep.join(EXAMPLES_GLOB))

path_args = []
            
for example_path in python_examples:
    if 'cooling_rate.py' in example_path:
        for i in range(4):
            path_args.append((example_path, i, None))
    elif 'freefall.py' in example_path:
        for i in range(2):
            path_args.append((example_path, None, i))
    else:
        path_args.append((example_path, None, None))

@pytest.mark.parametrize('example_path,primordial_chemistry,metal_cooling',
                         path_args)
def test_examples(example_path, primordial_chemistry, metal_cooling):
    env = dict(os.environ)
    if primordial_chemistry is not None:
        env['PRIMORDIAL_CHEM'] = str(primordial_chemistry)
    if metal_cooling is not None:
        env['METAL_COOLING'] = str(metal_cooling)
    python_executable = 'python'
    with temporary_directory() as tmpdir:
        command = '%s %s' % (python_executable, example_path)
        try:
            subprocess.check_output(
                command.split(' '), stderr=subprocess.STDOUT,
                cwd=tmpdir, env=env)
        except subprocess.CalledProcessError as er:
            raise RuntimeError('Command %s failed with return code %s '
                               'and the following output: %s' %
                               (command, er.returncode, er.output))

        example_base = os.path.split(example_path)[1].strip('.py')
        possible_file = '.'.join([example_base, 'h5'])

        if not os.path.exists(os.sep.join([tmpdir, possible_file])):
            return

        answer_path = os.sep.join([os.path.dirname(
            os.path.abspath(__file__)), 'example_answers'])
        if primordial_chemistry is not None:
            example_base = example_base + 'pc%s' % primordial_chemistry
        answer_name = '.'.join([example_base, 'h5'])

        ds_old = yt.load(os.sep.join([answer_path, answer_name]))
        ds_new = yt.load(os.sep.join([tmpdir, answer_name]))

        ad_old = ds_old.data
        ad_new = ds_new.data

        for field_name in ds_old.field_list:
            assert_allclose(ad_old[field_name].v, ad_new[field_name].v,
                            err_msg="Field mismatch: %s." % str(field_name))
