#!/usr/bin/env python3

########################################################################
#
# standalone python script that is used to help with testing that the
# fields of certain structs are dynamically accessible.
#
# This is NOT a direct part of the pytest framework
# - for compatability, we explicitly choose to:
#   - only use modules from the standard library
#   - only use features from python 3.7 or newer
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import argparse
import re
import sys
import warnings
import xml.etree.ElementTree  # may not be the optimal xml parser


def _find_element_by_ids(id_str, root, expect_single=False):
    # id_str provides a list of one or more space separated ids
    id_set = set(id_str.split(" "))
    assert len(id_set) or (len(id_set) != 0 and not expect_single)

    # implicit assumption: every element has a unique id
    out = [e for e in root if e.attrib["id"] in id_set]
    return out


def _field_type_props(elem, root):
    assert elem.tag == "Field"
    ptr_level = 0
    while True:
        elem = _find_element_by_ids(elem.attrib["type"], root, True)[0]
        if elem.tag == "PointerType":
            ptr_level += 1
        elif elem.tag == "CvQualifiedType":
            pass  # just ignore this
        else:
            return elem.attrib["name"] + (ptr_level * "*")


def _unwrap_struct_from_typedef(elem, root):
    # this is called when elem corresponds to a type declared with typedef
    # -> this function returns the xml-object representing the underlying struct
    #    that is aliased by the typedef
    _tmp = _find_element_by_ids(elem.attrib["type"], root, expect_single=True)[0]
    if _tmp.tag == "ElaboratedType":
        struct_elem = _find_element_by_ids(
            _tmp.attrib["type"], root, expect_single=True
        )[0]
    else:
        struct_elem = _tmp
    assert struct_elem.tag == "Struct"  # <-- sanity check!
    return struct_elem


def main(args):
    # castxml was previously used to parse a c++ file and spit out an xml file that
    # summarizes all of the contained information. We write the results to:
    out_f = args.output

    # this is the struct name we are hunting for:
    struct_name = args.type

    # parse the xml file
    tree = xml.etree.ElementTree.parse(args.input)

    root = tree.getroot()
    assert root.tag == "CastXML"
    m = re.match(
        r"^(?P<major>\d+)\.(?P<minor>\d+)\.(?P<patch>\d+)", root.attrib["format"]
    )
    try:
        version_info = (int(m["major"]), int(m["minor"]), int(m["patch"]))
    except Exception:
        raise RuntimeError(
            f"Unknown CastXML format version fmt: {root.attrib!r}"
        ) from None
    if (version_info < (1, 1, 0)) or (version_info > (1, 4, 0)):
        warnings.warn(
            "Code has only beeen tested against CastXML format versions 1.1.0 thru "
            f"1.4.0. The file produced by CastXML uses version {root.attrib['format']}"
        )

    # Since castxml is technically designed for C++ there will be some extra info that
    # we don't care about

    # under the root node, the children are listed in a flat structure can be
    # - typedef declarations, struct/class declarations, namespace declarations,
    #   Field declarations, class member declarations, global variable
    #   declarations
    # - it also has entries for describing other types

    # Now extract the necessary information
    # 1: find the element with name attribute that matches struct_name
    matches = tree.getroot().findall(
        ".//Struct[@name='" + struct_name + "']"
    ) + tree.getroot().findall(".//Typedef[@name='" + struct_name + "']")

    if len(matches) == 0:
        raise RuntimeError(f"no struct name '{struct_name}' was found")
    elif len(matches) > 2:
        raise RuntimeError(f"found more than 2 matches for '{struct_name}'")
    elif (len(matches) == 1) and (matches[0].tag == "Struct"):
        # in this case, we defined a struct named {struct_name} and didn't
        # use typedef to declare any type aliases of the struct
        # -> in other words, we can only refer to the type in our C code as
        #    `struct {struct_name}`
        struct_elem = matches[0]
    elif (len(matches) == 1) and (matches[0].tag == "Typedef"):
        # in this case, we used typedef to define a type called {struct_name}
        # that aliases an "anonymous struct"
        # -> in other words, we can only refer to the type in our C code as
        #    `{struct_name}`
        struct_elem = _unwrap_struct_from_typedef(matches[0], root)
    elif (
        (len(matches) == 2)
        and (matches[0].tag == "Struct")
        and (matches[1].tag == "Typedef")
    ):
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

    # 2: find the child elements of the struct
    out = []
    for member in _find_element_by_ids(
        struct_elem.attrib["members"], root, expect_single=False
    ):
        if member.tag == "Field":  # ignore auto-generated class methods
            data_member_name = member.attrib["name"]
            data_member_type = _field_type_props(member, root)
            out_f.write(f"{data_member_name}, {data_member_type}\n")
    return out


parser = argparse.ArgumentParser(
    prog="castxml_output_reader",
    description="makes a list of all members of a type from the output of castxml",
)
parser.add_argument("-i", "--input", required=True, help="path to xml file")
parser.add_argument("-t", "--type", required=True, help="name of the type")
parser.add_argument("-o", "--output", default=sys.stdout, help="optional output file")

if __name__ == "__main__":
    main(parser.parse_args())
