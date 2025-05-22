import os
import pytest
import sys
import yt

from numpy.testing import assert_allclose

from pygrackle.__config__ import _is_editable_installation
from pygrackle.utilities.model_tests import model_sets
from pygrackle.utilities.testing import run_command

from testing_common import grackle_python_dir

pytestmark = pytest.mark.skipif(
    not _is_editable_installation(),
    reason="this module currently requires an editable installation"
)

python_example_dir = os.path.join(grackle_python_dir, "examples")

# collect all of the python-examples and various model configurations into
# a list of tuples
all_sets = []
for model_name, model in model_sets.items():
    for par_index in range(len(model["parameter_sets"])):
        for input_index in range(len(model["input_sets"])):
            all_sets.append((model_name, par_index, input_index))

@pytest.mark.parametrize("model_name, par_index, input_index", all_sets)
def test_model(answertestspec, tmp_path, model_name, par_index, input_index):
    """
    Each execution tests a python example with a set of inputs

    Each of the parameters down below is a fixture

    Parameters
    ----------
    answertestspec: AnswerTestSpec
        A fixture that provides info about the answer-testing configuration
    tmp_path: pathlib.Path
        A custom built-in fixture provided by pytest that specifies a pre-made
        temporary directory that is named for the current test
    """

    if (model_name == "yt_grackle") and ("YT_DATA_DIR" not in os.environ):
        pytest.skip("YT_DATA_DIR env variable isn't defined")

    script_path = os.path.join(python_example_dir, f"{model_name}.py")
    command = f"{sys.executable} {script_path} {par_index} {input_index}"

    if True:
        rval = run_command(command, timeout=60, cwd=tmp_path)
        assert rval, "example didn't complete succesfully"

        output_basename = f"{model_name}_{par_index}_{input_index}.h5"
        output_file = os.path.join(str(tmp_path), output_basename)
        answer_path = os.path.join(answertestspec.answer_dir, output_basename)

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
