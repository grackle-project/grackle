########################################################################
#
# Test the API for dynamically accessing fields of chemistry_data
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import contextlib
import io
import os
import shutil
import sys

import numpy as np
import pytest


from pygrackle import setup_fluid_container, constants
from pygrackle.utilities.data_path import (
    _make_config_pair,
    _fnames_in_registry,
)
from pygrackle.utilities.grdata import main
from pygrackle.utilities.physical_constants import sec_per_Myr
from pygrackle.utilities.testing import assert_allequal_arraydict, ensure_dir

from test_query_units import _setup_generic_chemistry_data

# we probably don't have to skip everything
if not hasattr(os, "putenv"):
    pytest.skip(
        "several tests need os.putenv to work properly", allow_module_level=True
    )

# _ENV_VAR holds the list of environment variables that could affect the
# location of the data directory
if sys.platform.startswith("darwin"):
    _ENV_VARS = ("HOME", "GRACKLE_DATA_DIR")
else:
    _ENV_VARS = ("HOME", "GRACKLE_DATA_DIR", "XDG_DATA_HOME")


def _ensure_removed(d, key):
    try:
        del d[key]
    except KeyError:
        pass


@contextlib.contextmanager
def modified_env(new_env_vals, extra_cleared_variables=None):
    """
    Temporarily overwrite the environment variables. This is necessary to test C
    extensions that rely upon the environment variables
    """
    if extra_cleared_variables is None:
        extra_cleared_variables = None

    # record the original values for any variable we will overwrite
    original_vals = {}
    try:
        for var in filter(lambda e: e not in new_env_vals, extra_cleared_variables):
            original_vals[var] = os.environ.get(var, None)
            _ensure_removed(os.environ, var)

        for var, new_val in new_env_vals.items():
            original_vals[var] = os.environ.get(var, None)
            if new_val is None:
                _ensure_removed(os.environ, var)
            else:
                os.environ[var] = new_val

        yield

    finally:
        # restore to the initial values
        for var, val in original_vals.items():
            if val is None:
                _ensure_removed(os.environ, var)
            else:
                os.environ[var] = val


class DataFileManagementHarness:
    """
    This is a wrapper around the cli interface provided by pygrackle.

    This mainly exists to make it easier for us to wrap a standalone script
    in the future that isn't part of pygrackle
    """

    def __init__(self, config_pair=None):
        self.config_pair = config_pair
        self.fnames_in_registry = _fnames_in_registry()

    def __call__(self, args):
        """pass in cli args. The exit code and the captured stdout is returned"""
        if (args is None) or isinstance(args, str) or not isinstance(args[0], str):
            raise RuntimeError("invalid args sanity check failed!")
        config_pair = self.config_pair
        if config_pair is None:
            config_pair = _make_config_pair()
        tmp = io.StringIO()
        with contextlib.redirect_stdout(tmp):
            exitcode = main(*config_pair, prog_name="python -m pygrackle", args=args)
        return exitcode, tmp.getvalue().rstrip()

    def version_dir_path(self):
        rc, current_version_data_path = self(["getpath", "--vdata"])
        if rc != 0:
            raise RuntimeError("something went horribly wrong")
        return current_version_data_path

    def data_dir_path(self):
        rc, current_version_data_path = self(["getpath", "--data-dir"])
        if rc != 0:
            raise RuntimeError("something went horribly wrong")
        return current_version_data_path


# this one is locked environment variables as they are set right now
_static_GRDATA = DataFileManagementHarness(_make_config_pair())
# this one will be affected by changes in the environment variable
_flexible_GRDATA = DataFileManagementHarness()


@contextlib.contextmanager
def tmpversiondir_with_file(input_path, env, fname=None, *, cleanup_on_close=False):
    """
    A context manager that sets up a temporary on disk that appears (to
    the Grackle library), as if the grdata tool set up a data directory
    (the location is governed by the environment variables specified by
    env), that contains a single file called ``fname``, which is a copy
    of the has the file at `input_path`.

    In practice, the data-directory structure may not actually be
    managaed by the grdata tool. Consequently, some implementation
    details (e.g. related to deduplication) may not be defined. But,
    that's ok since the logic in the Grackle library should only care
    about whether a file (or link) shows up in the version directory.

    Parameters
    ----------
    input_path : str
        the path to the file we will copy
    env : dict of strs
        Dictionary holding the new values that we will use
    fname : Optional, str
        This is the name of the file as it appears inside of the
        versiondir. When not specified, this is inferred from
        input_path
    """
    for var, val in env.items():
        if var not in _ENV_VARS:
            raise ValueError(f"{var} isn't a known overridable env variable.")
    if not os.path.isfile(input_path):
        raise ValueError("input_path must specify a real file")

    with modified_env(env, extra_cleared_variables=_ENV_VARS):
        try:
            data_dir = _flexible_GRDATA.data_dir_path()
            if os.path.isdir(data_dir) and (len(os.listdir(data_dir)) > 0):
                raise ValueError(
                    "sanity check: this context manager requires that you specify "
                    "environment variables that lead to a data directory that doesn't "
                    "exist yet (or at least is empty)"
                )

            version_dir = _flexible_GRDATA.version_dir_path()
            if fname is None:
                # in the future, we may want to actually invoke GRDATA to make the copy
                # in this case
                fname = os.path.basename(input_path)

            ensure_dir(version_dir)
            full_path = os.path.join(version_dir, fname)
            shutil.copy(input_path, full_path)
            yield full_path

        finally:
            if cleanup_on_close:
                shutil.rmtree(_flexible_GRDATA.data_dir_path())


def _check_valid_datafile_fname(fname):
    if fname not in _static_GRDATA.fnames_in_registry:
        pytest.skip(
            f"test is broken since {fname} is not a datafile distributed "
            "with the current version Grackle"
        )


@pytest.fixture(scope="function")
def managed_datafile(request, tmp_path):
    """
    A pytest fixture that ensures that a data-directory (and associated
    environment variables specifying its location) is correctly configured
    so that Grackle's internal logic can automatically lookup the
    location of a standard datafile for the duration of a test.

    The standard datafile is "CloudyData_UVB=HM2012.h5".

    For the sake of convenience, this fixture passes provides the full
    path of the datafile (in that data directory) to the test.

    This operates in 2 modes:
      1. when `hasattr(request, "param", None)` is `None`, we use the
         existing data directory (essentially we ignore the tmp_path
         fixture).
      2. otherwise, we use the `tmpversiondir_with_file` context manager
         to temporarily (for the duration of the test) delete any/all
         environment variables that could control the location of the
         data-directory and replace it with the environment variable
         specified by `param.request`.
         - That environment variable hints at the location of a
           temporary data directory.
         - The location of that directory is controlled by the path
           provided by the pytest's `tmp_path` fixture.
         - We also copy the standard datafile into the appropriate
           location within the data directory so that the test can
           actually read in the data file.

    Note
    ----
    If we want to parameterize the actual name of the file, then maybe we should return
    some kind of factory?
    """

    fname = "CloudyData_UVB=HM2012.h5"
    _check_valid_datafile_fname(fname)

    existing_fname_path = os.path.join(_static_GRDATA.version_dir_path(), fname)

    if getattr(request, "param", None) is None:
        full_path = existing_fname_path
        yield full_path
    else:
        env_var = request.param
        with tmpversiondir_with_file(
            input_path=existing_fname_path, env={env_var: str(tmp_path)}
        ) as full_path:
            yield full_path

def setup_generic_problem(parameter_overrides={}):
    """set up a really simplistic problem"""
    chem = _setup_generic_chemistry_data(
        initial_redshift=2.7, parameter_overrides=parameter_overrides
    )
    # the precise details don't really matter here...
    dt = sec_per_Myr / chem.time_units
    fc = setup_fluid_container(
        chem,
        density=1.67e-24,
        temperature=np.geomspace(1e3, 1e7, num=11),
        metal_mass_fraction=0.01,  # kinda arbitrary
        state="ionized",
        converge=False,
    )
    return fc, dt


@pytest.mark.parametrize(
    "managed_datafile",
    ([pytest.param(None, id = "default-datadir")] +
     [pytest.param(var, id=f"arbitrary-{var}") for var in _ENV_VARS]),
    indirect=True
)
def test_autofile_equivalence(managed_datafile):
    """
    A parameterized test that confirms that grackle produces the same
    exact result (for a generic test problem) when you:
    - you pass grackle_data_file a full path to the data file
    - automatic lookup is used to infer the full path (to the same file)

    This test uses a parameterized fixture that may
    - use the existing data directory variable,
    - or use a custom environment variable that specifies the location
      of the data file (in this case, the variable points to a location in
      a temporary directory, where the datafile has been copied to)

    Essentially, the use of parametrized fixtures let us confirm that
    Grackle's internal logic searches for the data files in the right
    locations.
    """

    full_path = managed_datafile
    fname = os.path.basename(full_path)

    assert os.path.isfile(full_path)  # sanity check

    # generate a simple test problem
    fc_ref, dt = setup_generic_problem(
        parameter_overrides={"grackle_data_file": full_path}
    )
    fc_ref.solve_chemistry(dt)

    # rerun the same problem, but now don't use the full path
    fc_other, _ = setup_generic_problem(
        parameter_overrides={
            "grackle_data_file": fname,
            "grackle_data_file_options": constants.GR_DFOPT_MANAGED
        }
    )
    fc_other.solve_chemistry(dt)
    assert_allequal_arraydict(fc_ref, fc_other)


def test_autofile_fail_unknown_file():
    # verify that the autofile machinery properly tells Grackle to abort initialization
    # when we specify an invalid filename
    chem = _setup_generic_chemistry_data(
        initial_redshift=0.0,
        skip_initialize=True,
        parameter_overrides={
            "grackle_data_file": "not-a-file.png",
            "grackle_data_file_options": constants.GR_DFOPT_MANAGED
        },
    )
    assert chem.initialize() == constants.GR_FAIL


def test_autofile_fail_known_missing_file(tmp_path):
    # verify that the autofile machinery properly tells Grackle to abort initialization
    # when we specify a filename known to Grackle but that is missing

    fname_to_copy = "CloudyData_UVB=HM2012.h5"
    alt_fname = "CloudyData_UVB=FG2011.h5"
    _check_valid_datafile_fname(fname_to_copy)
    _check_valid_datafile_fname(alt_fname)

    file_to_copy = os.path.join(_static_GRDATA.version_dir_path(), fname_to_copy)

    with tmpversiondir_with_file(
        input_path=file_to_copy,
        env={"GRACKLE_DATA_DIR": str(tmp_path)},
    ):
        chem = _setup_generic_chemistry_data(
            initial_redshift=0.0,
            skip_initialize=True,
            parameter_overrides={
                "grackle_data_file": alt_fname,
                "grackle_data_file_options": constants.GR_DFOPT_MANAGED
            },
        )
        assert chem.initialize() == constants.GR_FAIL


def test_autofile_fail_bad_checksum(tmp_path):
    # verify that the autofile machinery properly tells Grackle to abort initialization
    # when we specify a filename known to Grackle, that exists, but has the wrong
    # checksum value

    fname_to_copy = "CloudyData_UVB=HM2012.h5"
    alt_fname = "CloudyData_UVB=FG2011.h5"
    _check_valid_datafile_fname(fname_to_copy)
    _check_valid_datafile_fname(alt_fname)

    file_to_copy = os.path.join(_static_GRDATA.version_dir_path(), fname_to_copy)

    # for this test, we intentionally copy a file and give it the wrong name
    with tmpversiondir_with_file(
        input_path=file_to_copy,
        env={"GRACKLE_DATA_DIR": str(tmp_path)},
        fname=alt_fname,
    ):
        chem = _setup_generic_chemistry_data(
            initial_redshift=0.0,
            skip_initialize=True,
            parameter_overrides={
                "grackle_data_file": alt_fname,
                "grackle_data_file_options": constants.GR_DFOPT_MANAGED
            },
        )
        assert chem.initialize() == constants.GR_FAIL

