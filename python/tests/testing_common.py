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

from gracklepy.utilities.misc import dirname

grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
grackle_data_dir = os.path.join(grackle_install_dir, "input")
grackle_python_dir = os.path.join(grackle_install_dir, "src", "python")
