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
import warnings
import xml.etree.ElementTree # may not be the optimal xml parser

import pytest

from pygrackle.grackle_wrapper import _wrapped_c_chemistry_data

_CASTXML_INSTALLED = shutil.which('castxml') is not None

def _find_element_by_ids(id_str, root, expect_single = False):
    # id_str provides a list of one or more space separated ids
    id_set = set(id_str.split(' '))
    assert len(id_set) or (len(id_set) != 0 and not expect_single)

    # implicit assumption: every element has a unique id
    out = [e for e in root if e.attrib['id'] in id_set]
    return out

def _field_type_props(elem, root):
    assert elem.tag == 'Field'
    ptr_level = 0
    while True:
        elem = _find_element_by_ids(elem.attrib['type'], root, True)[0]
        if elem.tag == 'PointerType':
            ptr_level += 1
        elif elem.tag == 'CvQualifiedType':
            pass # just ignore this
        else:
            return elem.attrib['name'] + (ptr_level * "*")

def _unwrap_struct_from_typedef(elem, root):
    # this is called when elem corresponds to a type declared with typedef
    # -> this function returns the xml-object representing the underlying struct
    #    that is aliased by the typedef
    _tmp = _find_element_by_ids(elem.attrib['type'], root,
                                expect_single = True)[0]
    if _tmp.tag == 'ElaboratedType':
        struct_elem = _find_element_by_ids(_tmp.attrib['type'], root,
                                           expect_single = True)[0]
    else:
        struct_elem = _tmp
    assert struct_elem.tag == 'Struct'  # <-- sanity check!
    return struct_elem

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
        # tell castxml about include-paths
        '-I', os.path.dirname(path),
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

    if not os.path.isfile(xml_fname):
        raise RuntimeError(
            f"something went wrong while executing {' '.join(command)}"
        )

    # Step 4: parse the output and then delete the temporary file
    tree = xml.etree.ElementTree.parse(xml_fname)
    os.remove(xml_fname)

    # in version 1.1.0 of gccxml format, there is a root node whose tag is
    # CastXML & lists the version attribute
    root = tree.getroot()
    assert root.tag == 'CastXML'
    if root.attrib['format'] not in ['1.1.0', '1.1.5']:
        warnings.warn(
            "Code was only tested against CastXML format versions 1.1.0 & "
            "1.1.5. The file produced by CastXML uses version "
            f"{root.attrib['format']}"
        )

    # under the root node, the children are listed in a flat structure can be
    # - typedef declarations, struct/class declarations, namespace declarations,
    #   Field declarations, class member declarations, global variable
    #   declarations
    # - it also has entries for describing other types

    # Step 5: now extract the necessary information
    # Step 5a: find the element with name attribute that matches struct_name
    matches = (
        tree.getroot().findall(".//Struct[@name='" + struct_name + "']") +
        tree.getroot().findall(".//Typedef[@name='" + struct_name + "']")
    )

    if len(matches) == 0:
        raise RuntimeError(f"no struct name '{struct_name}' was found")
    elif len(matches) > 2:
        raise RuntimeError(f"found more than 2 matches for '{struct_name}'")
    elif (len(matches) == 1) and (matches[0].tag == 'Struct'):
        # in this case, we defined a struct named {struct_name} and didn't
        # use typedef to declare any type aliases of the struct
        # -> in other words, we can only refer to the type in our C code as
        #    `struct {struct_name}`
        struct_elem = matches[0]
    elif (len(matches) == 1) and (matches[0].tag == 'Typedef'):
        # in this case, we used typedef to define a type called {struct_name}
        # that aliases an "anonymous struct"
        # -> in other words, we can only refer to the type in our C code as
        #    `{struct_name}`
        struct_elem = _unwrap_struct_from_typedef(matches[0], root)
    elif ((len(matches) == 2) and (matches[0].tag == 'Struct') and
          (matches[1].tag == 'Typedef')):
        # in this case we probably defined a struct called {struct_name} and
        # used typedef to declare an alias type called {struct_name}
        # -> in other words, the types `struct {struct_name}` and
        #    `{struct_name}` refer to the same type
        struct_elem = matches[0]
        # perform a sanity check to confirm that typedef indeed aliases the
        # corresponding struct
        # -> note: I'm not entirely sure this is always guaranteed to work
        assert struct_elem == _unwrap_struct_from_typedef(matches[1], root)
    else:
        raise RuntimeError("SOMETHING WENT WRONG")

    # Step 5b, find the child elements of the struct
    out = []
    for member in _find_element_by_ids(struct_elem.attrib['members'], root,
                                       expect_single = False):
        if member.tag == 'Field': # ignore auto-generated class methods
            out.append(
                (_field_type_props(member, root), member.attrib['name'])
            )
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
    for dtype, name in member_list:
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
