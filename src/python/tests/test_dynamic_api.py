########################################################################
#
# Test the API for dynamically accessing fields of chemistry_data
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from numpy.testing import assert_raises

from pygrackle import chemistry_data

from pygrackle.grackle_wrapper import _wrapped_c_chemistry_data

def _test_dynamic_api(param_list, ok_vals, bad_vals, expected_output_type):
    # test that the chemistry_data() class supports access to an attribute
    # named for element in key_list.
    #
    # most other tests of chemistry_data() also implicitly test the dynamic_api

    def _check_expected_type(param, val):
        if not isinstance(val, expected_output_type):
            raise AssertionError(
                f"expected '{param}' to be an instance of "
                f"'{expected_output_type.__name__}', not an instance of "
                f"'{type(val).__name__}'"
            )

    obj = chemistry_data()
    for param in param_list:
        val = getattr(obj, param)
        _check_expected_type(param, val)

        for val_to_set in ok_vals:
            setattr(obj, param, val_to_set)
            val = getattr(obj, param)
            _check_expected_type(param, val)

            # check equality between val_to_set and val
            if isinstance(val_to_set, str) and expected_output_type is bytes:
                assert val.decode('ascii') == val_to_set
            else:
                assert val == val_to_set

        # confirm that setting param to invalid values will raise TypeError
        for val_to_set in bad_vals:
            with assert_raises(TypeError):
                setattr(obj, param, val_to_set)

def test_dynamic_api_int():
    _test_dynamic_api(_wrapped_c_chemistry_data.int_keys(),
                      [342], [2.5, 'my string'], int)

def test_dynamic_api_double():
    _test_dynamic_api(_wrapped_c_chemistry_data.double_keys(),
                      [342, 2.0], ['my string'], float)

def test_dynamic_api_string():
    _test_dynamic_api(_wrapped_c_chemistry_data.string_keys(),
                      [b'dummy_bytes', 'dummy_str'], [1, 3.0], bytes)
