import importlib
import json
import numpy as np
from matplotlib import pyplot as plt
import os
import pytest
import shutil
import sys
import yt

from numpy.testing import assert_allclose

from gracklepy.__config__ import _is_editable_installation
from gracklepy.utilities.model_tests import model_parametrization

from testing_common import grackle_python_dir

pytestmark = pytest.mark.skipif(
    not _is_editable_installation(),
    reason="this module currently requires an editable installation"
)

def _get_model_main_fn(module_name):
    # the implementation is inspired by
    # https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    file_path = os.path.abspath(
        os.path.join(grackle_python_dir, "examples", f"{module_name}.py"))
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module.main


model_fns = {
    model : _get_model_main_fn(module_name=model)
    for model in ["cooling_cell", "cooling_rate", "freefall", "yt_grackle"]
}

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

def _infer_symlog_linthresh(y):
    # infer kwarg based on the y-values being plotted for ax.set_yscale("symlog", ...)
    # -> By default, we try to pick values that minimize the size of the region
    #    where the y-axis crosses 0.
    # -> a value of None means that we give up and just want builtin defaults

    if np.ndim(y) != 1:
        return None
    abs_y = np.abs(y)
    tmp = abs_y[(abs_y != 0) & np.isfinite(abs_y)]
    if tmp.size == 0:
        return None
    linthresh = tmp.min()
    return linthresh if linthresh < np.nanmax(abs_y) else None

def _subplot_text(ax, msg):
    ax.text(
        0.5,
        0.5,
        msg,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes
    )

def _careful_plot_logy(ax, x, y, **kw):
    if "linestyle" in kw:
        raise NotImplementedError("not currently supported")
    out = []

    has_positive = (y > 0).any()
    if has_positive:
        out.extend(ax.plot(x, y, linestyle = "-", **kw))

    if (y < 0).any():
        _kw = kw
        if has_positive and "label" in kw:
            _kw = {k: v for k, v in kw.items() if k != "label"}
        out.extend(ax.plot(x, -1*y, linestyle = "--", **_kw))
    ax.set_yscale("log")
    return out


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
        json.dump(my_pars, f, indent=2)

    output_lines = ["\nMax relative differences:"]
    do_diff = ds1.data["data", my_ivar].size == ds2.data["data", my_ivar].size

    # prepare properties of each output figure
    _height_ratios = [3, 1, 1]
    _height_inches = 4.8 * sum(_height_ratios) / _height_ratios[0]
    subplots_kw = {
        "nrows": len(_height_ratios),
        "ncols": 1,
        "sharex": True,
        "figsize": (6.4, _height_inches),
        "gridspec_kw": {"height_ratios": _height_ratios}
    }

    for field in ds1.field_list:
        if field[1] in _exclude_fields:
            continue
        if field[1] == my_ivar:
            continue

        if (ds1.data[field] == 0).all() and (ds2.data[field] == 0).all():
            continue

        fig, (ax1, ax2, ax3) = plt.subplots(**subplots_kw)

        new_kw = {"color": "red", "label": "new", "alpha": 0.8}
        old_kw = {"color": "green", "label": "old", "alpha": 0.8}

        if field[1] == "cooling_rate":
            lt1 = _infer_symlog_linthresh(ds1.data[field].ndview)
            lt2 = _infer_symlog_linthresh(ds2.data[field].ndview)
            if lt1 is None and lt2 is None:
                ax1.set_yscale("symlog")  # <- just use defaults!
            else:
                lt = min(
                    np.inf if lt1 is None else lt1, np.inf if lt2 is None else lt2
                )
                ax1.set_yscale("symlog", linthresh=lt)
            ax1.plot(ds1.data["data", my_ivar], ds1.data[field], **new_kw)
            ax1.plot(ds2.data["data", my_ivar], ds2.data[field], **old_kw)
        else:
            _careful_plot_logy(
                ax1, ds1.data["data", my_ivar], ds1.data[field], **new_kw
            )
            _careful_plot_logy(
                ax1, ds2.data["data", my_ivar], ds2.data[field], **old_kw
            )
        xunits = getattr(ds1.field_info["data", my_ivar], "units", "")
        yunits = getattr(ds1.field_info[field], "units", "")
        ax3.set_xlabel(f"{my_ivar} [{xunits}]")
        ax1.set_ylabel(f"{field[1]} [{yunits}]")
        ax2.set_ylabel("relative difference")
        ax3.set_ylabel("abs difference")
        ax1.legend(loc="best")
        ax3.set_xscale("log")

        if do_diff:
            abs_diff = np.abs(ds1.data[field] - ds2.data[field])
            diff = abs_diff / ds2.data[field]
            my_line = f"\t{field[1]}: {diff.max().d}."
            output_lines.append(my_line)

            for (ax, cur_diff) in [(ax2, diff), (ax3, abs_diff)]:
                if (cur_diff > 0).any():
                    ax.loglog(ds1.data["data", my_ivar], cur_diff, color="black")
                elif (cur_diff == 0).all():
                    _subplot_text(ax, "All elements are identical")
                else:
                    _subplot_text(ax, "All differences are zeros or NaNs")

        else:

            for ax in [ax2,ax3]:
                _subplot_text(
                    ax, "Unable to compare: data arrays have different sizes"
                )

        fig.tight_layout()
        fig.savefig(os.path.join(test_dir, f"{field[1]}.png"))
        plt.close(fig)

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

    return_val = model_fns[model_name](
        args=[
            f"--out-dir={tmp_path!s}",
            model_variant,
            str(par_index),
            str(input_index)
        ],
    )

    model_par = f"{model_name}_{model_variant}_{par_index}_{input_index}"
    assert return_val == 0, f"Model {model_par} didn't complete succesfully."

    output_basename = f"{model_par}.h5"
    output_file = os.path.join(str(tmp_path), output_basename)
    answer_path = os.path.join(answertestspec.answer_dir, output_basename)

    if answertestspec.generate_answers:
        # use shutil.move over os.rename in case output_file & answer_path
        # specify locations on different file systems
        # -> this is common since we using pytest's tmp_path machinery to
        #    determine output_file. Unless overriden, this use's the
        #    platform's standard directory for temporary files
        # -> On some platforms, especially linux, this location is entirely
        #    stored in RAM
        shutil.move(output_file, answer_path)
    else:
        assert os.path.exists(answer_path), \
            f"No gold-standard answer file exists at {answer_path}"

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
