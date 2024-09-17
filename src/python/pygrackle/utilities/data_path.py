########################################################################
#
# Try to get a path to the Grackle input data.
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os

from pygrackle.__config__ import _is_editable_installation
from pygrackle.utilities.misc import dirname

grackle_data_dir = os.environ.get("GRACKLE_DATA_DIR")
if (grackle_data_dir is None) and _is_editable_installation():
    # Note, this only works with an editable install of pygrackle.
    _install_dir = dirname(os.path.abspath(__file__), level=5)
    grackle_data_dir = os.path.join(_install_dir, "input")
elif (grackle_data_dir is None):
    raise RuntimeError(
        "in non-editable pygrackle installations, like this one, "
        f"grackle_data_dir can only be imported from {__file__} if it is set "
        "by the GRACKLE_DATA_DIR environment variable"
    )

if not os.path.isdir(grackle_data_dir):
    raise RuntimeError(
        f"grackle_data_dir not set to a valid directory: {grackle_data_dir}. "
        "Use the GRACKLE_DATA_DIR environment variable to set path to "
        "Grackle data.")
