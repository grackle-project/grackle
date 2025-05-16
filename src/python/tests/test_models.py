import numpy as np
from matplotlib import pyplot as plt
import os
import pytest
import sys
import yaml
import yt

from numpy.testing import assert_allclose

from pygrackle.utilities.model_tests import model_parametrization
from pygrackle.utilities.testing import run_command

from testing_common import grackle_python_dir

python_example_dir = os.path.join(grackle_python_dir, "examples")

_ivars = {
    "cooling_rate": "temperature",
    "freefall": "density",
    "cooling_cell": "time",
}

def get_ivar(my_test):
    for key, ivar in _ivars.items():
        if key in my_test:
            return ivar
    return None

_exclude_fields = [f"{ax}_velocity" for ax in "xyz"]

def compare_model_results(compare_dir, model_par, ds1, ds2):
    """
    Do additional comparison between models.

    This will create a folder containing plots of fields
    comparing the existing answers and the ones just generated.
    """

    my_ivar = get_ivar(model_par)
    if my_ivar is None:
        return

    test_dir = os.path.join(compare_dir, model_par)
    os.makedirs(test_dir, exist_ok=True)

    notes_fn = os.path.join(test_dir, "notes.txt")
    with open(notes_fn, mode="w") as f:
        my_pars = {k: str(v) for k, v in ds1.parameters.items()}
        yaml.dump(my_pars, stream=f)

    output_lines = ["\nMax relative differences:"]
    do_diff = ds1.data["data", my_ivar].size == ds2.data["data", my_ivar].size

    for field in ds1.field_list:
        if field[1] in _exclude_fields:
            continue
        if field[1] == my_ivar:
            continue

        if (ds1.data[field] == 0).all() and (ds2.data[field] == 0).all():
            continue

        has_positive = False
        if (ds1.data[field] > 0).any():
            plt.loglog(ds1.data["data", my_ivar], ds1.data[field],
                       color="red", label="new", alpha=0.8)
            plt.loglog(ds2.data["data", my_ivar], ds2.data[field],
                       color="green", label="old", alpha=0.8)
            has_positive = True

        if (ds1.data[field] < 0).any():
            if has_positive:
                label1 = label2 = None
            else:
                label1 = "new"
                label2 = "old"
            plt.loglog(ds1.data["data", my_ivar], -ds1.data[field],
                       color="red", linestyle="--", label=label1, alpha=0.8)
            plt.loglog(ds2.data["data", my_ivar], -ds2.data[field],
                       color="green", linestyle="--", label=label2, alpha=0.8)

        xunits = getattr(ds1.field_info["data", my_ivar], "units", "")
        yunits = getattr(ds1.field_info[field], "units", "")
        plt.xlabel(f"{my_ivar} [{xunits}]")
        plt.ylabel(f"{field[1]} [{yunits}]")
        plt.legend(loc="best")

        if do_diff:
            diff = np.abs(ds1.data[field] - ds2.data[field]) / ds2.data[field]
            my_line = f"\t{field[1]}: {diff.max().d}."
            output_lines.append(my_line)

            if (diff > 0).any():
                plt.twinx()
                plt.loglog(ds1.data["data", my_ivar], diff,
                           color="black", alpha=0.8)
                plt.ylabel("relative difference")

        plt.tight_layout()
        plt.savefig(os.path.join(test_dir, f"{field[1]}.png"))
        plt.clf()

    if not do_diff:
        output_lines = ["\nData arrays have different sizes."]

    output_lines.append("")
    with open(notes_fn, mode="a") as f:
        f.write("\n".join(output_lines))

@pytest.mark.parametrize("model_name, model_variant, par_index, input_index",
                         model_parametrization)
def test_model(answertestspec, tmp_path, model_name, model_variant,
               par_index, input_index):
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
    command = f"{sys.executable} {script_path} {model_variant} {par_index} {input_index}"

    model_par = f"{model_name}_{model_variant}_{par_index}_{input_index}"
    rval = run_command(command, timeout=60, cwd=tmp_path)
    assert rval, f"Model {model_par} didn't complete succesfully."

    output_basename = f"{model_par}.h5"
    output_file = os.path.join(str(tmp_path), output_basename)
    answer_path = os.path.join(answertestspec.answer_dir, output_basename)

    if answertestspec.generate_answers:
        os.rename(output_file, answer_path)
    else:
        assert os.path.exists(answer_path)

        ds1 = yt.load(output_file)
        ds2 = yt.load(answer_path)
        assert str(ds1.parameters["format_version"]) == \
          str(ds2.parameters["format_version"])
        assert ds1.field_list == ds2.field_list

        if answertestspec.model_comparison_dir is not None:
            compare_model_results(answertestspec.model_comparison_dir,
                                  model_par, ds1, ds2)

        for field in ds1.field_list:
            err_msg = f"Model mismatch: {model_name}, {model_variant}, " + \
              f"{par_index}, {input_index}: {field}."
            assert_allclose(ds1.data[field], ds2.data[field],
                            atol=0, rtol=1e-8, err_msg=err_msg)
