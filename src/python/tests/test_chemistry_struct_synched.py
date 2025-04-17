########################################################################
#
# Explicitly test that the API for dynamically accessing fields of
# chemistry_data is synchronized with the members of chemistry_data
# (this is meant to identify the scenario where a new member gets added to
# chemistry_data but the dynamic API is not synchronized)
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os
import os.path
import shutil
import subprocess
import tempfile

import pytest

from pygrackle.grackle_wrapper import _wrapped_c_chemistry_data
from testing_common import grackle_install_dir

_CASTXML_INSTALLED = shutil.which('castxml') is not None

def query_struct_fields(struct_name, path):
    """
    Query the members of a struct using the castxml commandline tool.

    castxml is used to parse a file c++ (specified by path) and spit out an xml
    file that summarizes all of the contained information. Since it's
    technically designed for C++ there will be some extra info that we don't
    care about

    Note
    ----
    pygccxml is a python module that can be used to simplify this function.
    However, it doesn't look like that module has been updated recently. To try
    to mitigate future potential module compatibility issues, we currently
    choose not to use it.
    """
    assert _CASTXML_INSTALLED

    # Step 1: get temp file name (it's probably ok to we use old insecure API)
    xml_fname = tempfile.mktemp(suffix='.xml')

    # Step 2: build up the command
    command = [
        'castxml',
        # tell castxml's compiler to preprocess/compile but don't link
        '-c',
        # tell castxml's compiler that the language is c++
        '-x', 'c++',
        # the next required option configure castxml's internal Clang compiler.
        # The second part needs to specify an installed compiler (we may need
        # to support configuration of that path)
        '--castxml-cc-gnu', 'g++',
        # specify the output xml format version
        '--castxml-output=1',
        # specify the output location
        '-o', xml_fname,
        # specify the src file to process
        path
    ]

    # Step 3: run the command
    subprocess.run(command, check=True)

    # Step 4: run script to analyze the outputs
    command = [
        f"{grackle_install_dir}/tests/scripts/castxml_output_reader.py",
        "--input", xml_fname,
        "--type", struct_name
    ]
    rslt = subprocess.run(command, capture_output=True, check=True, encoding="utf8")
    out = []
    for line in rslt.stdout.splitlines():
        if len(line) == "":
            continue
        pairs = line.split(",")
        assert len(pairs) == 2
        name, dtype = pairs
        out.append((name.strip(), dtype.strip()))
    return out

@pytest.mark.skipif(not _CASTXML_INSTALLED,
                    reason = 'requires the castxml program')
def test_grackle_chemistry_field_synched():
    # use the castxml to construct a list of the members of the chemistry_data
    # struct (by directly parsing the header file)
    member_list = query_struct_fields(
        struct_name = "chemistry_data",
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "../../include/grackle_chemistry_data.h")
    )

    # now, categorize the fields by their datatype
    field_sets = {'char*' : set(), 'int' : set(), 'double' : set()}
    for name, dtype in member_list:
        field_sets[dtype].add(name)

    # finally lets compare
    for param_type, ref_set in [('int', field_sets['int']),
                                ('double', field_sets['double']),
                                ('string', field_sets['char*'])]:
        if param_type == 'int':
            parameters = _wrapped_c_chemistry_data.int_keys()
        elif param_type == 'double':
            parameters = _wrapped_c_chemistry_data.double_keys()
        elif param_type == 'string':
            parameters = _wrapped_c_chemistry_data.string_keys()
        else:
            raise RuntimeError(f"unrecognized parameter type: {param_type}")

        diff = ref_set.symmetric_difference(parameters)
        for parameter in diff:
            if parameter in ref_set:
                raise RuntimeError(
                    f"{parameter}, a {param_type} field of the chemistry_data "
                    "struct is not accessible from the dynamic api"
                )
            else:
                raise RuntimeError(
                    f"the dynamic api provides access to '{parameter}', a "
                    f"{param_type} parameter. But it's not a member of the "
                    "chemistry_data struct."
                )
