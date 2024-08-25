########################################################################
#
# miscellaneous utilities
#
#
# Copyright (c) Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from numpy import VisibleDeprecationWarning
import os
import warnings

def dirname(path, level=1):
    """
    Multi-level version of os.path.dirname.
    """
    if not isinstance(level, int) or level < 1:
        raise ValueError(
            f"level must be a positive integer: {level}.")
    for i in range(level):
        path = os.path.dirname(path)
    return path

def issue_deprecation_warning(msg, stacklevel=3):
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)
