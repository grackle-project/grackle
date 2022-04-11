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

from numpy import \
    VisibleDeprecationWarning

import warnings

def issue_deprecation_warning(msg, stacklevel=3):
    warnings.warn(msg, VisibleDeprecationWarning, stacklevel=stacklevel)
