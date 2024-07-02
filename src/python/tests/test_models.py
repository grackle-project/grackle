import os
import pytest
import sys

from pygrackle.utilities.model_tests import model_sets
from pygrackle.utilities.testing import \
    dirname, \
    run_command, \
    temporary_directory

grackle_install_dir = dirname(os.path.abspath(__file__), level=4)
python_example_dir = os.path.join(
    grackle_install_dir, "src", "python", "examples")

generate_results = \
  int(os.environ.get("GENERATE_MODEL_TEST_RESULTS", 0)) == 1

all_sets = []
for model_name, model in model_sets.items():
    for par_index in range(len(model["parameter_sets"])):
        for input_index in range(len(model["input_sets"])):
            all_sets.append((model_name, par_index, input_index))

@pytest.mark.parametrize("model_name, parameter_index, input_index", all_sets)
def test_model(model_name, parameter_index, input_index):
    script_path = os.path.join(python_example_dir, f"{model_name}.py")
    command = f"{sys.executable} {script_path} {parameter_index} {input_index}"
    with temporary_directory():
        rval = run_command(command, timeout=60)
        assert rval
