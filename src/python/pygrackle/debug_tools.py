########################################################################
#
# Provides tools to aide with debugging
#
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

# we probably don't want to import anything from here into the main __init__.py
# file in order to avoid introducing the a hard dependence on h5py

import contextlib

import h5py

from .fluid_container import (
    FluidContainer,
    _convert_field_data_member_name_to_key
)
from .utilities.init_from_dicts import (
    _configure_chemistry_data,
    _configure_codeunits
)

_META_ATTR_PREFIX = "@:?GRACKLEDUMP_"
# keep the following 2 synchronized with the C source code
_META_ATTRS = {"protocol_version" : "@:?GRACKLEDUMP_VERSION",
               "complete_dump" : "@:?GRACKLEDUMP_DONE"}
assert all(a.startswith(_META_ATTR_PREFIX) for a in _META_ATTRS.values())

def _load_h5dump_grp(grp, omit_empty_dset = False):
    data = {}
    protocol_version = None
    complete_dump = None
    for key in grp.attrs:
        if key == _META_ATTRS["protocol_version"]:
            protocol_version = grp.attrs[key]
        elif key == _META_ATTRS["complete_dump"]:
            complete_dump = grp.attrs[key]
        elif key.startswith(_META_ATTR_PREFIX):
            raise RuntimeError(f"unrecognized metadata: {key!r}")
        else:
            data[key] = grp.attrs[key]

    for key in grp:
        if key.startswith(_META_ATTR_PREFIX):
            raise RuntimeError(f"The {key!r} group/dataset is named like "
                               "a metadata attribute")
        if isinstance(grp[key], h5py.Group):
            data[key] = _load_h5dump_grp(grp[key], omit_empty_dset)
        elif grp[key].shape is None: # in this case, the dset is empty
            if not omit_empty_dset:
                data[key] = None
        else:
            data[key] = grp[key][...]

    if protocol_version is None:
        raise RuntimeError("there is no protocol version")
    elif protocol_version != 1:
        raise RuntimeError(
            "we support protocol version 1, not {protocol_version}")
    elif complete_dump is None:
        raise RuntimeError("the datadump is incomplete")
    elif complete_dump != 1:
        raise RuntimeError("something went wrong when the datadump was "
                           "originally written")
    return data

def load_h5dump_dicts(file, *, omit_empty_dset = False):
    """
    Returns the relevant Grackle-information

    Parameters
    ----------
    file : path-like or h5py.Group
        Either a path to an HDF5 file or an already openned HDF5 file (or group
        within that file). This file should have been created with grackle's
        debugging C API.
    omit_empty_dset : bool, optional
        By default, `None` is used to represent any empty datasets (which
        correspond to NULL pointer arrays). When True, they are skipped
        alltogether.

    Notes
    -----
    The layout of the HDF5 file format is NOT stable and will change over time.
    """
    if isinstance(file, h5py.Group): # note: h5py.File subclasses h5py.File
        cm = contextlib.nullcontext(file)
    else:
        cm = h5py.File(file, 'r')

    with cm as f:
        if "grackle_statedump" not in f:
            raise ValueError("The specified file must have a grackle_statedump "
                             "group.")
        out = _load_h5dump_grp(f["grackle_statedump"], omit_empty_dset)

        # sanity-check
        assert all(isinstance(e, dict) for e in out.values())
    return out

def _load_chemistry_data_from_h5data(data, *,
                                     override_grackle_data_file = None,
                                     cur_units_usage = 'default'):
    # helper function for load_chemistry_data_from_h5dump

    # for conciseness, abbreviate the standard unit keys
    _INIT_U_KEY, _CUR_U_KEY = "initial_code_units", "current_code_units"

    if cur_units_usage == 'default':
        init_u, overwrite_u = _INIT_U_KEY, _CUR_U_KEY
    elif cur_units_usage == 'exclusive':
        init_u, overwrite_u = _INIT_U_KEY, _CUR_U_KEY
    elif cur_units_usage == 'exclusive-fallback':
        init_u = _INIT_U_KEY if (_INIT_U_KEY in data) else _CUR_U_KEY
        overwrite_u = _CUR_U_KEY
    else:
        raise ValueError("the cur_units_usage kwarg must be one of "
                         "'default', 'exclusive', 'exclusive-fallback'")

    # now we check if the needed data was actually saved
    for key in ("chemistry_data", init_u, overwrite_u):
        if key not in data:
            raise RuntimeError(f"{key} struct was not dumped to disk")

    # initialize Grackle
    # -> we require that the inputs exhaustively specify all attributes as a
    #    sort of sanity-check (it's a cheap extra check to ensure that might
    #    catch changes in the c-library that were not communicated to this
    #    machinery)
    chem_data_dict = data["chemistry_data"]
    if override_grackle_data_file is not None:
        chem_data_dict = data["chemistry_data"].copy()
        chem_data_dict['grackle_data_file'] = override_grackle_data_file
    my_chem = _configure_chemistry_data(chem_data_dict, data[init_u],
                                        require_exhaustive = True)

    # overwrite the units stored by Grackle (if necessary)
    # -> again, we require exhaustive attributes as a sanity check!
    if init_u != overwrite_u:
        _configure_codeunits(units_obj = my_chem,
                             code_units_config = data[overwrite_u],
                             require_exhaustive = True)
    return my_chem

def load_chemistry_data_from_h5dump(f, *, override_grackle_data_file = None,
                                    cur_units_usage = 'default'):
    """
    Loads the fully initialized chemistry_data instance from an h5dump

    Parameters
    ----------
    f : path-like or h5py.Group
        Either a path to an HDF5 file or an already openned HDF5 file (or group
        within that file). This file should have been created with grackle's
        debugging C API.
    override_grackle_data_file : str, optional
        When not `None`, this parameter is used to override the
        grackle_data_file chemistry parameter.
    cur_units_usage : str, optional
        Policy for specifying how the loaded current-code-units are used.
        Should be one of:

            - `'default'`: This is the default behavior and is intended to
              reproduce the existing behavior from the code executed before the
              dump. Grackle with initial-code-units and then overwrite
              with current-code-units before returning. This is
            - `'exclusive'`: The current code-units are exclusively used (for
              initialization and there is no overwriting). In other words,
              the initial code-units are never used.
            - `'exclusive-fallback'`: The `'default'` behavior is adopted, if
              possible. Otherwise, the `'exclusive'` behavior is adopted.

    Notes
    -----
    The layout of the HDF5 file format is NOT stable and will change over time.
    """
    data = load_h5dump_dicts(f)
    return _load_chemistry_data_from_h5data(
        data, override_grackle_data_file = override_grackle_data_file,
        cur_units_usage = cur_units_usage)

# I actually don't want to make this a part of the API... but, enough
# transformations need to be made on the data to actually warrant making this
# part of the API

def load_FluidContainer_from_h5dump(f, *, override_grackle_data_file = None,
                                    cur_units_usage = 'default',
                                    flatten_field_data = False):
    """
    Loads the fully initialized FluidCollection instance from an h5dump

    Parameters
    ----------
    f : path-like or h5py.Group
        Either a path to an HDF5 file or an already openned HDF5 file (or group
        within that file). This file should have been created with grackle's
        debugging C API.
    override_grackle_data_file : str, optional
        When not `None`, this parameter is used to override the
        grackle_data_file chemistry parameter.
    cur_units_usage : str, optional
        Policy for specifying how the loaded current-code-units are used.
        Should be one of:

            - `'default'`: This is the default behavior and is intended to
              reproduce the existing behavior from the code executed before the
              dump. Grackle with initial-code-units and then overwrite
              with current-code-units before returning. This is
            - `'exclusive'`: The current code-units are exclusively used (for
              initialization and there is no overwriting). In other words,
              the initial code-units are never used.
            - `'exclusive-fallback'`: The `'default'` behavior is adopted, if
              possible. Otherwise, the `'exclusive'` behavior is adopted.
    flatten_field_data : bool, optional
        When True, we will flatten the shapes of the grackle field data before
        using it.

    Notes
    -----
    The layout of the HDF5 file format is NOT stable and will change over time.
    """
    data = load_h5dump_dicts(f, omit_empty_dset = False)
    my_chem = _load_chemistry_data_from_h5data(
        data, override_grackle_data_file = override_grackle_data_file,
        cur_units_usage = cur_units_usage)

    _KEY = "grackle_field_data"
    if _KEY not in data:
        raise RuntimeError("grackle_field_data struct wasn't dumped to disk")

    field_data = data[_KEY]
    fc = FluidContainer(my_chem, n_vals = field_data["density"].size)
    for member, vals in field_data.items():
        if member in ['grid_dx', 'grid_rank', 'grid_dimension', 'grid_end',
                      'grid_start']:
            continue

        key = _convert_field_data_member_name_to_key(member)

        if vals is None:
            # member corresponds to an empty-dataset. As a consistency-check
            # (to help end-users identify issues), it would be nice to confirm
            # whether or not the field is needed by the current calculation.
            # Unfortunately, FluidContainer adopts a strategy where it may
            # define unnecessary arrays.
            continue
        elif key not in fc:
            raise RuntimeError("SOMETHING IS WRONG: FluidContainer doesn't "
                               f"have a {key!r} key")
        fc[key][...] = vals
    return fc

