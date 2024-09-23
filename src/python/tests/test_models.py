import os
import pytest
import sys
import yt

from numpy.testing import assert_allclose

from pygrackle.utilities.model_tests import model_sets
from pygrackle.utilities.testing import \
    run_command, \
    temporary_directory

from testing_common import grackle_python_dir

python_example_dir = os.path.join(grackle_python_dir, "examples")

all_sets = []
for model_name, model in model_sets.items():
    for par_index in range(len(model["parameter_sets"])):
        for input_index in range(len(model["input_sets"])):
            all_sets.append((model_name, par_index, input_index))

@pytest.mark.parametrize("model_name, par_index, input_index", all_sets)
def test_model(answertestspec, model_name, par_index, input_index):
    script_path = os.path.join(python_example_dir, f"{model_name}.py")
    command = f"{sys.executable} {script_path} {par_index} {input_index}"
    with temporary_directory():
        rval = run_command(command, timeout=60)
        assert rval

        output_file = f"{model_name}_{par_index}_{input_index}.h5"
        answer_path = os.path.join(answertestspec.answer_dir, output_file)

        if answertestspec.generate_answers:
            os.rename(output_file, answer_path)
        else:
            assert os.path.exists(answer_path)

            ds1 = yt.load(output_file)
            ds2 = yt.load(answer_path)
            ds1.parameters["format_version"] == ds2.parameters["format_version"]
            assert ds1.field_list == ds2.field_list

            for field in ds1.field_list:
                err_msg = f"Model mismatch: {model_name}, {par_index}, " + \
                  f"{input_index}: {field}."
                assert_allclose(ds1.data[field], ds2.data[field],
                                atol=0, rtol=1e-8, err_msg=err_msg)
