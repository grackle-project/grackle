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

from pygrackle.utilities.misc import dirname

generate_test_results = \
  int(os.environ.get("GENERATE_PYGRACKLE_TEST_RESULTS", 0)) == 1

# Note, this only works for an editable pygrackle installation,
# but the assumption is that one will only be running the tests
# from an editable install.
grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
grackle_data_dir = os.path.join(grackle_install_dir, "input")
grackle_python_dir = os.path.join(grackle_install_dir, "src", "python")
test_answers_dir = os.path.join(grackle_python_dir, "tests", "test_answers")
