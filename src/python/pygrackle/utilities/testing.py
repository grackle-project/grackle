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

def _valid_single_cmp_kwargs(d):
    return all(k in {'rtol', 'atol', 'equal_nan'} and not isinstance(v, dict)
               for k, v in d.items())

def _with_default_cmp_kwargs(kwargs):
    tmp = {'rtol' : 0.0, 'atol' : 0.0, 'equal_nan' : False}
    tmp.update(kwargs)
    return tmp

def assert_allclose_resultdict(actual, reference, cmp_kwarg_spec = dict(),
                               err_msg=''):
    """
    Raises an AssertionError if two objects are not equal up to desired
    tolerance.

    Parameters
    ----------
    actual : mapping
         A mapping of arrays obtained in a calculation
    reference : mapping
         A mapping of reference arrays
    cmp_kwarg_spec : dict, optional
         This is used to specify the comparison properties for arrays in actual
         and reference. By default, this requires exact equivalence. You can
         also use this to specify a single set of kwargs for `np.isclose` (i.e.
         `'rtol'`, `'atol'`, `'equal_nan'`) to be used 
         for comparing ALL arrays. Alternatively, you can use this to a dict
         of kwargs (e.g. a dict of dicts). In this last case, kwargs associated
         a given key value, `k`, will be used as kwargs for the comparison
         between `actual[k]` and `reference[k]`. In this last case, if `actual`
         has a key not in `cmp_kwarg_spec`, the function will try to make use
         of any kwargs associated with the `None` key.
    err_msg : str
         Custom error message to be printed in case of failure.
    """
    keys = _fetch_keys(actual, reference, err_msg = err_msg)

    # define the _get_kw function to retrieve 
    if not isinstance(cmp_kwarg_spec, dict):
        raise TypeError("cmp_kwarg_spec must be None or a dict")
    elif _valid_single_cmp_kwargs(cmp_kwarg_spec):
        _get_cmp_kw = lambda key: _with_default_cmp_kwargs(cmp_kwarg_spec)
    else:
        for k,v in cmp_kwarg_spec.items(): # error check
            if not _valid_single_cmp_kwargs(v):
                raise ValueError(
                    f"the entry associated with cmp_kwarg_spec[{k!r}] doesn't "
                    "specify a valid comparison kwarg set for isclose")
            elif (k is not None) and (k not in keys):
                raise ValueError(
                    "cmp_kwarg_spec specifies a comparison kwargs for a key, "
                    f"{k!r} that doesn't exist in actual or reference")
        def _get_cmp_kw(key):
            if key in cmp_kwarg_spec:
                return _with_default_cmp_kwargs(cmp_kwarg_spec[key])
            return _with_default_cmp_kwargs(cmp_kwarg_spec.get(None, {}))

    for key in keys:
        # it would be nice to customize this quite a bit more. We could:
        # -> optionally customize the error-message header (we could let caller
        #    specify a callback to provide more details about what exactly the
        #    data held by key represents)
        # -> optionally customize reporting of the actual error itself
        #    (particularly if the array is relatively large). It would be nice
        #    to report the index and values of the arrays at:
        #      - the first location of nonclose vals
        #      - the location of nonclose vals where the magnitude of rtol is
        #        maximized
        #      - the location of nonclose vals where the magnitude of rtol is
        #        maximized
        cmp_kwargs = _get_cmp_kw(key)
        np.testing.assert_allclose(
            actual = actual[key], desired = reference[key],
            err_msg = (
                f"(In comparison of arrays associated with the {key!r} key)\n"
                f"{err_msg}"),
            verbose = True, **cmp_kwargs)

