########################################################################
#
# Common variables for testing.
#
#
# Copyright (c), Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################


import os

from pygrackle.__config__ import _is_editable_installation
from pygrackle.utilities.misc import dirname

if _is_editable_installation():
    grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
    grackle_data_dir = os.path.join(grackle_install_dir, "input")
    grackle_python_dir = os.path.join(grackle_install_dir, "src", "python")
else:
    raise RuntimeError(
        "the current version of pygrackle is not an editable installation. No "
        "tests that must import from {__file__} can be run"
    )
