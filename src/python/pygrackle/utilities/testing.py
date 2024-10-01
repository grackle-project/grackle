########################################################################
#
# Testing imports
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import contextlib
import importlib
import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal, \
    assert_approx_equal, assert_array_almost_equal, assert_equal, \
    assert_array_less, assert_string_equal, assert_array_almost_equal_nulp,\
    assert_allclose, assert_raises
import os
import shutil
import subprocess
import tempfile

def assert_rel_equal(a1, a2, decimals, err_msg='', verbose=True):
    if isinstance(a1, np.ndarray):
        assert(a1.size == a2.size)
        # Mask out NaNs
        a1[np.isnan(a1)] = 1.0
        a2[np.isnan(a2)] = 1.0
    elif np.any(np.isnan(a1)) and np.any(np.isnan(a2)):
        return True
    return assert_almost_equal(np.array(a1)/np.array(a2), 1.0,
                               decimals, err_msg=err_msg,
                               verbose=verbose)

def random_logscale(log_min, log_max, size=1, random_state=None):
    if random_state is None:
        random_state = np.random.default_rng()
    log_val = (log_max - log_min) * random_state.random(size) + log_min
    return np.power(10, log_val)

def requires_module(module):
    """
    Decorator that takes a module name as an argument and tries to import it.
    If the module imports without issue, the function is returned, but if not,
    a null function is returned. This is so tests that depend on certain modules
    being imported will not fail if the module is not installed on the testing
    platform.
    """
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        return func
    try:
        importlib.import_module(module)
    except ImportError:
        return ffalse
    else:
        return ftrue

def run_command(command, timeout=None):
    try:
        proc = subprocess.run(command, shell=True, timeout=timeout)
        if proc.returncode == 0:
            success = True
        else:
            success = False
    except subprocess.TimeoutExpired:
        print ("Process reached timeout of %d s. (%s)" % (timeout, command))
        success = False
    except KeyboardInterrupt:
        print ("Killed by keyboard interrupt!")
        success = False
    return success

@contextlib.contextmanager
def temporary_directory():
    curdir = os.getcwd()
    tmpdir = tempfile.mkdtemp(dir=curdir)
    os.chdir(tmpdir)
    try:
        yield tmpdir
    finally:
        os.chdir(curdir)
        shutil.rmtree(tmpdir)

def ensure_dir(path):
    r"""Parallel safe directory maker."""
    if os.path.exists(path):
        return path

    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
        else:
            raise
    return path
