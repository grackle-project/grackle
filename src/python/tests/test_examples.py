
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
import shutil
import subprocess
import tempfile
import yt

from matplotlib.testing.compare import \
    compare_images
from pygrackle.utilities.testing import \
    assert_equal, \
    assert_almost_equal

current_path = os.path.abspath(__file__)

EXAMPLES_GLOB = [os.path.dirname(
    os.path.dirname(current_path)), 'examples', '*.py']

python_examples = glob.glob(os.sep.join(EXAMPLES_GLOB))


@contextlib.contextmanager
def temporary_directory():
    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=curdir)
    try:
        yield tmpdir
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)


def example_test(example_path, primordial_chemistry=None):
    env = dict(os.environ)
    if primordial_chemistry is not None:
        env['PRIMORDIAL_CHEM'] = str(primordial_chemistry)
    if 'PYTHON_EXE' in env:
        python_executable = env['PYTHON_EXE']
    else:
        python_executable = 'python'
    with temporary_directory() as tmpdir:
        output = ''
        try:
            print tmpdir
            output = subprocess.check_output(
                [python_executable, example_path], stderr=subprocess.STDOUT,
                cwd=tmpdir, env=env)
        except subprocess.CalledProcessError as er:
            command = '%s %s' % (python_executable, example_path)
            raise RuntimeError('Command %s failed with return code %s '
                               'and the following output: %s' %
                               (command, er.returncode, output))

        example_base = os.path.split(example_path)[1].strip('.py')
        possible_files = [
            '.'.join([example_base, ext]) for ext in ['h5', 'png']]

        if not all([os.path.exists(os.sep.join([tmpdir, pf]))
                    for pf in possible_files]):
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
            assert_almost_equal(ad_old[field_name], ad_new[field_name])

        answer_name = '.'.join([example_base, 'png'])

        result = compare_images(os.sep.join([answer_path, answer_name]),
                                os.sep.join([tmpdir, answer_name]), 12)
        assert_equal(result, None, result)


def test_examples():
    for example_path in python_examples:
        if 'cooling_rate.py' in example_path:
            for i in range(4):
                yield example_test, example_path, i
        else:
            yield example_test, example_path
