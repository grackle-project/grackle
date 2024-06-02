########################################################################
#
# routines to help initialize the chemistry_data extension type from
# dictionaries
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from pygrackle import chemistry_data
from pygrackle.grackle_wrapper import (
    _wrapped_c_chemistry_data,
    _CODE_UNITS_ATTR_SET
)

# The functions defined in this file should not be used outside of
# Grackle. (since they may change at any time).

def _copy_dict_to_attrs(dest, src_dict, short_attr_descr, *,
                        expected_attrs = None, require_exhaustive = False):
    """
    Helper function that uses the entries of a dict to assign values to
    attributes of a generic object

    Parameters
    ----------
    dest
        The destination object whose attributes will be modified.
    src_dict : dict
        The dict containing the key-value pairs that will be assigned.
    short_attr_descr : str
        A short generic description of any given attribute in dest. If errors
        arise this is used as a singular noun in the error message.
    expected_attrs : container, optional
        specifies the set of keys that src_dict is allowed to initialize
    require_exhaustive : bool, optional
        When True, an error is raised when `src_dict` doesn't contain a key
        for every value in expected_attrs.
    """
    if (expected_attrs is None) and require_exhaustive:
        raise ValueError("require_exhaustive can only be true when "
                         "expected_attrs is not None")
    elif expected_attrs is None:
        try:
            for k,v in src_dict.items():
                setattr(dest, k, v)
        except AttributeError:
            raise ValueError(
                f"unknown {short_attr_descr} was specified: {k!r}") from None
    else:
        for k,v in src_dict.items():
            if k not in expected_attrs:
                raise ValueError(
                    f"Unknown {short_attr_descr} was specified: {k!r}")
            setattr(dest, k, v)

        if require_exhaustive and (len(src_dict) > len(expected_attrs)):
            for k in expected_attrs:
                if k not in src_dict:
                    raise ValueError(f"missing {short_attr_descr}: {k!r}")

def _configure_codeunits(units_obj, code_units_config, *,
                         require_exhaustive = False):
    """
    Setup values related to the code_units struct.
    """
    _copy_dict_to_attrs(dest = units_obj, src_dict = code_units_config,
                        short_attr_descr = "code-units-attribute",
                        expected_attrs = _CODE_UNITS_ATTR_SET,
                        require_exhaustive = require_exhaustive)
    units_obj.set_velocity_units()


_CHEM_PARAMETERS = frozenset(_wrapped_c_chemistry_data.int_keys() +
                             _wrapped_c_chemistry_data.double_keys() +
                             _wrapped_c_chemistry_data.string_keys())

def _configure_chemistry_data(chemistry_data_parameters, code_units, *,
                              require_exhaustive = False):
    """
    Construct a fully initialized `pygrackle.grackle_wrapper.chemistry_data`
    instance from the 2 dictionaries

    Parameters
    ----------
    chemistry_data_parameters : dict
        Specifies the parameter names held by Grackle's `chemistry_data` struct
        and their associated values
    code_units : dict
        Specifies the attributes held by Grackle's `code_units` struct and their
        associated values.
    require_exhaustive : bool, optional
        When True, an error is raised when `chemistry_data_parameters` and
        `code_units` don't exhaustively specify all values held by Grackle's
        `chemistry_data` and `code_units` structs.
    """

    # this is pretty much the control flow in a c or c++ program

    # step 1: declare/initialize chemistry_data and code_units structs (they
    #         are both wrapped by `pygrackle.grackle_wrapper.chemistry_data`)
    my_chemistry = chemistry_data()

    # step 2: configure parameters that are stored in chemistry_data
    _copy_dict_to_attrs(dest = my_chemistry,
                        src_dict = chemistry_data_parameters,
                        short_attr_descr = "chemistry-data-parameter",
                        expected_attrs = _CHEM_PARAMETERS,
                        require_exhaustive = require_exhaustive)

    # step 3: configure the stuff related to the code_units struct
    _configure_codeunits(units_obj = my_chemistry,
                         code_units_config = code_units,
                         require_exhaustive = require_exhaustive)

    # step 4: actually initialize the underlying chemistry_data struct
    if my_chemistry.initialize() != 1:
        raise RuntimeError("There was an issue with initialization")

    return my_chemistry

