#!/usr/bin/env python3
"""
This forwards the arguments on to castxml, if the castxml is installed. Otherwise
it puts an empty output file at the specified location
"""

import shutil
import subprocess
import sys

_PROG = "castxml"


def _get_output_location(args):
    for i, arg in enumerate(args):
        if arg == "-o" and i + 1 < len(args):
            return args[i + 1]


if __name__ == "__main__":
    if shutil.which(_PROG) is None:
        output_location = _get_output_location(sys.argv)
        if output_location is None:
            sys.exit(1)
        f = open(output_location, "w")  # <- truncates file if it exists
        f.close()
        sys.exit(0)
    else:
        code = subprocess.run([_PROG] + sys.argv[1:]).returncode
        sys.exit(code)
