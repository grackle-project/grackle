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

from gracklepy.__config__ import _is_editable_installation
from gracklepy.utilities.misc import dirname

if _is_editable_installation():
    # Note, this only works with an editable install of gracklepy.
    _install_dir = dirname(os.path.abspath(__file__), level=5)
    grackle_data_dir = os.path.join(_install_dir, "input")
else:
    raise RuntimeError(
        "in non-editable gracklepy installations, like this one, "
        f"grackle_data_dir cannot be imported from {__file__}."
    )
