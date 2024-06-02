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

import importlib
import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal, \
    assert_approx_equal, assert_array_almost_equal, assert_equal, \
    assert_array_less, assert_string_equal, assert_array_almost_equal_nulp,\
    assert_allclose, assert_raises

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
        random_state = np.random.RandomState()
    log_val = (log_max - log_min) * random_state.random_sample(size) + log_min
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


def _fetch_keys(actual, reference, err_msg = ""):
    # check consistency in dictionary keys
    refkeys = reference.keys()
    refkey_set = set(refkeys)
    mismatch_keys = refkey_set.symmetric_difference(actual.keys())

    if len(mismatch_keys):
        shared_keys = list(refkey_set.intersection(actual.keys()))
        extra_ref, extra_actual = [], []
        for k in mismatch_keys:
            if k in refkeys:
                extra_ref.append(k)
            else:
                extra_actual.append(k)

        raise AssertionError(
            "The results are not equal to specified tolerance.\n"
            f"{err_msg}\n"
            "There is a keys mismatch. Both results have the keys:\n"
            f" {shared_keys!r}\n"
            "Extra Keys:\n"
            f" actual:    {extra_actual}\n"
            f" reference: {extra_ref}"
        )
    return list(refkeys)

def assert_allequal_arraydict(actual, reference, err_msg=''):
    """
    Raises an AssertionError if any contents of the 2 compared mappings of
    arrays are not EXACTLY equal

    Parameters
    ----------
    actual : mapping
         A mapping of arrays obtained in a calculation
    reference : mapping
         A mapping of reference arrays
    err_msg : str
         Custom error message to be printed in case of failure.

    Note
    ----
    A separate function is proposed as part of PR #195 to do approximate
    equality checks (like np.testing.assert_allclose).
    """
    keys = _fetch_keys(actual, reference, err_msg = err_msg)
    for key in keys:
        assert_array_equal(actual[key], reference[key], err_msg = err_msg,
                           strict = True)

