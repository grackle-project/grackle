########################################################################
#
# Test the command line tool managing test files
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
import hashlib
import io
import operator
import os
import shutil
import subprocess
import sys
from textwrap import indent
from typing import Any, NamedTuple

import pytest

# a goal here is to be able to run this test without installing pygrackle!
# -> in the near future, we will install grdata as a standalone command-line
#    script and it would really nice to be able to test the command-line script
#    without installing pygrackle
# -> when that time comes, we will modify the logic within the cli_app fixture
# -> currently, we need to include the following import. But, in the future, we
#    could add a new subcommand to grdata to make it unnecessary
from pygrackle.utilities.grdata import _parse_file_registry


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


@contextlib.contextmanager
def custom_datadir(path):
    """
    A contextmanager used to put the data directory at an arbitrary location.
    """
    clear_env = [var for var in _ENV_VARS if var != "GRACKLE_DATA_DIR"]
    with modified_env({"GRACKLE_DATA_DIR": path}, clear_env):
        yield


class GRDataExecErr(Exception):
    pass


class CLIApp:
    """exists to wrap the command-line-interface.

    We use this so that we can eventually test the application when it is
    configured as a standalone script.
    """

    def __init__(self, common_args, *, use_function_call=None):
        self.common_args = common_args
        self.use_function_call = use_function_call

    def __call__(self, subcommand_args, *, expect_success=None):
        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)

        all_args = self.common_args + subcommand_args
        if self.use_function_call is None:
            tmp = subprocess.run(all_args, capture_output=True)
            returncode = tmp.returncode
            stdout = tmp.stdout.decode("ascii")
            stderr = tmp.stderr.decode("ascii")
        else:
            fn = self.use_function_call
            with contextlib.ExitStack() as stack:
                f_out, f_err = [
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO())),
                    stack.enter_context(contextlib.redirect_stderr(io.StringIO())),
                ]
                try:
                    returncode = fn(all_args)
                except SystemExit as err:
                    returncode = err.code
                stdout = f_out.getvalue()
                stderr = f_err.getvalue()

        expected_result = (
            (expect_success is None)
            or (returncode == 0 and expect_success)
            or (returncode != 0 and not expect_success)
        )
        if not expected_result:
            detail_indent = "   >"
            msg_lines = [
                "Invocation of grdata produced an unexpected result:\n",
                f"  expected: {('failure', 'success')[expect_success]}\n",
                f"  args: {all_args}\n",
                "  env:\n",
            ]
            for var in _ENV_VARS:
                msg_lines.append(
                    f"{detail_indent}{var!r}: {os.environ.get(var,'<unset>')!r}\n"
                )
            msg_lines.append(f"  returncode: {tmp.returncode}\n")

            for stream, val in [("stdout", stdout), ("stderr", stderr)]:
                if val is None or len(val) == 0:
                    msg_lines.append(f"  {stream}: <N/A>\n")
                else:
                    msg_lines += [
                        f"  {stream}:\n",
                        indent(val.decode("ascii"), detail_indent),
                    ]
            raise GRDataExecErr("".join(msg_lines))
        return returncode, stdout.rstrip()

    def fetch(
        self,
        src_dir=None,
        *,
        file_list=None,
        untracked_dest_dir=None,
        expect_success=None,
    ):
        # a number of tests care about whether the command fails. (In some case we
        # actually expect it to fail). We return whether or not it was succesful.

        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)

        if src_dir is None:
            subcommand_args = ["fetch"]
        else:
            subcommand_args = ["fetch", "--from-dir", src_dir]

        if file_list is not None:
            subcommand_args += file_list

        if untracked_dest_dir is not None:
            subcommand_args += ["--untracked-dest-dir", untracked_dest_dir]
        return self(subcommand_args, expect_success=expect_success)[0] == 0

    def rm_vdata(self, version, *, omit_force=False, expect_success=None):
        # remove a whole version-directory. returns whether this was successful

        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)

        if version is None:
            # delete the associated version
            subcommand_args = ["rm", "--force", "--vdata"]
        elif isinstance(version, str):
            subcommand_args = ["rm", "--force", "--vdata", version]
        else:
            # this particular mistake occurs a surprising amount
            raise TypeError("version must be None or a str")
        if omit_force:
            subcommand_args.remove("--force")

        return self(subcommand_args, expect_success=expect_success)[0] == 0

    def rm_datastore(self, *, omit_force=False, expect_success=None):
        # remove a whole data-store. returns whether this was successful

        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        subcommand_args = ["rm", "--force", "--data-store"]
        if omit_force:
            subcommand_args.remove("--force")
        return self(subcommand_args, expect_success=expect_success)[0] == 0

    def showknownreg(self):
        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        return self(["showknownreg"], expect_success=True)[1]

    def calcreg(self, cksum_alg, dir_path):
        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        return self(
            ["calcreg", "--hash-name", cksum_alg, dir_path], expect_success=True
        )[1]

    def cksum_alg(self):
        # have pytest hide certain kinds of noisy tracebacks
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        return self(["--cksum-alg"], expect_success=True)[1]

    def version_dir_path(self):
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        return self(["getpath", "--vdata"], expect_success=True)[1]

    def data_dir_path(self):
        __tracebackhide__ = operator.methodcaller("errisinstance", GRDataExecErr)
        return self(["getpath", "--data-dir"], expect_success=True)[1]


@pytest.fixture(scope="module")
def cli_app():
    if False:
        assert sys.executable is not None
        return CLIApp([sys.executable, "-m", "pygrackle"])
    else:
        import pygrackle.__main__

        return CLIApp([], use_function_call=pygrackle.__main__.main)


_SHASUM_INSTALLED = shutil.which("shasum") is not None


def _calc_ref_cksum(contents, cksum_alg):
    if _SHASUM_INSTALLED:
        _algs = {"sha1": "1", "sha256": "256"}
        args = ["shasum", "--algorithm", _algs[cksum_alg], "-"]
        rslt_str = (
            subprocess.run(
                args, input=contents.encode("ascii"), check=True, capture_output=True
            )
            .stdout.rstrip()
            .decode("utf8")
        )
        if rslt_str.endswith("  -"):
            cksum = rslt_str[:-3].lower()
        else:
            raise RuntimeError(f"the output of shasum was unexpected: '{rslt_str}'")
    else:
        if isinstance(contents, str):
            contents = contents.encode("ascii")
        hash_obj = hashlib.new(cksum_alg)
        hash_obj.update(contents)
        cksum = hash_obj.hexdigest()
    return f"{cksum_alg}:{cksum}"


class DummyFileSpec(NamedTuple):
    contents_str: str
    sha1: str
    sha256: str


def _dummy_file_contents(variant=1, trailing_content=None):
    assert variant >= 0 and int(variant) == variant
    newline_str = "\n" * (variant + 1)
    contents_str = f"I am a test-file.{newline_str}Variant number {variant}\n"
    if trailing_content is not None:
        contents_str = contents_str + trailing_content
    return DummyFileSpec(
        contents_str,
        _calc_ref_cksum(contents_str, "sha1"),
        _calc_ref_cksum(contents_str, "sha256"),
    )


# here we define file-sets. Each fileset is a sequence of pairs specifying a filename
# and its contents. The idea is to act like these correspond to different grackle
# versions and make sure we can handle them appropriately

_DUMMY_SET_PRIMARY = (
    ("file-0.txt", _dummy_file_contents(1)),
    ("file-1.txt", _dummy_file_contents(2)),
)

_DUMMY_SET_RENAME = (
    _DUMMY_SET_PRIMARY[0],
    ("renamed-file-2.txt", _DUMMY_SET_PRIMARY[1][1]),
)

# this scenario shouldn't come up in practice (we replaced a file with a different one,
# of the same name), but we should still handle it properly
_DUMMY_SET_REPLACE = (
    _DUMMY_SET_PRIMARY[0],
    ("file-2.txt", _dummy_file_contents(2, trailing_content="version 2 of file\n")),
)


class FileSetTuple(NamedTuple):
    """Holds an object for each of the filesets"""

    # corresponds to the primary file-set
    primary: Any
    # exactly like the primary file-set, but the 2nd file was renamed
    rename: Any
    # exactly like the primary file-set, but the contents of the second file was changed
    replace: Any

    def get(self, key):
        if key in self._fields:
            return getattr(self, key)
        raise KeyError(key)


_DUMMY_SET_TUPLE = FileSetTuple(
    _DUMMY_SET_PRIMARY, _DUMMY_SET_RENAME, _DUMMY_SET_REPLACE
)


class DummyFileRepository(NamedTuple):
    test_dir: str  # the path reserved for the user to do stuff in

    # following variables specify properties for primary set of dummy files
    registry_path: FileSetTuple
    src_file_dir: FileSetTuple

    def cli_app_with_overrides(self, ref, kind, version_override=None):
        new_args = ["--testing-override-registry-file", self.registry_path.get(kind)]

        if version_override is not None:
            if not isinstance(version_override, str):
                raise TypeError("version_override must be a str")
            new_args += ["--testing-override-version-grackle", version_override]
        return CLIApp(
            ref.common_args + new_args, use_function_call=ref.use_function_call
        )


@pytest.fixture
def dummy_file_repo(tmp_path, cli_app):
    test_dir = os.path.join(tmp_path, "test-dir")
    os.mkdir(test_dir)

    path = os.path.join(tmp_path, "fixture_dir")

    cksum_kind = cli_app.cksum_alg()

    registry_paths, src_file_dirs = {}, {}

    for kind in _DUMMY_SET_TUPLE._fields:
        file_set = _DUMMY_SET_TUPLE.get(kind)

        registry_path = os.path.join(path, f"{kind}_file_registry.txt")
        src_file_dir = os.path.join(path, f"{kind}-ref-files")
        os.makedirs(src_file_dir)

        src_file_dirs[kind] = src_file_dir
        registry_paths[kind] = registry_path

        pairs = []
        for i, (fname, file_spec) in enumerate(file_set):
            full_path = os.path.join(src_file_dir, fname)
            with open(full_path, "w") as f:
                f.write(file_spec.contents_str)
            pairs.append((fname, getattr(file_spec, cksum_kind)))
        with open(registry_path, "w") as f:
            print(
                *["".join(['{"', p[0], '", "', p[1], '"}']) for p in pairs],
                sep=",\n",
                file=f,
            )

    yield DummyFileRepository(
        test_dir=test_dir,
        registry_path=FileSetTuple(**registry_paths),
        src_file_dir=FileSetTuple(**src_file_dirs),
    )


def test_showknownreg(dummy_file_repo, cli_app):
    # essentially, we are checking that the testing override works
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary")
    full_registry_str = app.showknownreg().rstrip()

    with open(dummy_file_repo.registry_path.primary, "r") as f:
        ref_full_registry_str = f.read().rstrip()
    assert full_registry_str == ref_full_registry_str


def test_calcreg(dummy_file_repo, cli_app):
    for alg in ["sha1", "sha256"]:
        registry_str = cli_app.calcreg(alg, dummy_file_repo.src_file_dir.primary)
        registry = _parse_file_registry(io.StringIO(registry_str))

        for i, (fname, cksum) in enumerate(sorted(registry.items())):
            if not hasattr(_DUMMY_SET_PRIMARY[i][1], alg):
                raise RuntimeError("This should never happen. Unclear what went wrong")
            ref = getattr(_DUMMY_SET_PRIMARY[i][1], alg)
            if cksum != ref:
                raise AssertionError(
                    f"calculation of the {alg} checksum for the dummy-file, {fname}, "
                    "may have revealed an issue in the command line tool's "
                    "internal checksum logic\n"
                    f"expected: {ref}\nactual: {cksum}"
                )


def _get_lockfile_path(datadir_path):
    return os.path.join(datadir_path, "lockfile")


def _get_datastore_dir(datadir_path):
    return os.path.join(datadir_path, "data-store-v1")


def _get_version_dir(datadir_path, version):
    return os.path.join(datadir_path, "data-store-v1", version)


def _get_managed_file(datadir_path, version, fname):
    return os.path.join(datadir_path, "data-store-v1", version, fname)


def _dummy_errmsg_writer(msg):
    return AssertionError(msg)


def check_version_data_dir_contents(
    version_dir_path,
    file_set,
    *,
    exhaustive_file_set=True,
    errmsg_writer=_dummy_errmsg_writer,
):
    __tracebackhide__ = True  # suppress noisy pytest tracebacks
    for fname, file_spec in file_set:
        full_path = os.path.join(version_dir_path, fname)

        if not os.path.isfile(full_path):
            raise errmsg_writer(
                f"file, {full_path}, doesn't exist after the last invocation of grdata"
            )
        with open(full_path, "r") as f:
            contents = f.read()
        if contents != file_spec.contents_str:
            raise errmsg_writer(
                f"the file, {full_path}, doesn't have the correct contents"
            )

    if exhaustive_file_set and (len(os.listdir(version_dir_path)) != len(file_set)):
        raise errmsg_writer(
            f"the directory, {version_dir_path}, doesn't contain the right number of "
            "entries."
        )


def check_version_data_dir(
    datadir_path, version, lockfile_should_exist=False, file_set=None, *, err_msg=None
):
    __tracebackhide__ = True  # suppress noisy pytest tracebacks
    required_paths = [
        ("data directory", datadir_path),
        ("user-data directory", _get_version_dir(datadir_path, version)),
        ("version data directory", os.path.join(datadir_path, "user-data")),
    ]

    def prep_assertion_err(nominal_err):
        if err_msg is None:
            return AssertionError(nominal_err)
        return AssertionError(f"{err_msg}\n\n{nominal_err}")

    for descr, path in required_paths:
        if not os.path.isdir(path):
            raise prep_assertion_err(
                f"the {descr}, {path} does not exist following an invocation of grdata"
            )

    lockfile_path = _get_lockfile_path(datadir_path)
    if lockfile_should_exist and not os.path.isfile(lockfile_path):
        raise prep_assertion_err(
            f"the lockfile should exist at {lockfile_path} following "
            "the invocation of grdata."
        )
    elif (not lockfile_should_exist) and os.path.isfile(lockfile_path):
        raise prep_assertion_err(
            "a lockfile shouldn't exist after the last invocation of grdata."
        )

    exhaustive_file_set = True
    if file_set is None:
        file_set = []
        exhaustive_file_set = False
    check_version_data_dir_contents(
        version_dir_path=_get_version_dir(datadir_path, version),
        file_set=file_set,
        exhaustive_file_set=exhaustive_file_set,
        errmsg_writer=prep_assertion_err,
    )


def _check_removal(data_dir, version, retains_datastore):
    version_dir = _get_version_dir(data_dir, version)
    datastore_dir = _get_datastore_dir(data_dir)
    if os.path.isdir(version_dir):
        raise AssertionError(
            f"after a successful remove operation, the version-dir, {version_dir}, "
            "shouldn't exist"
        )
    elif os.path.isdir(datastore_dir) != retains_datastore:
        raise AssertionError(
            f"the data-store directory, {datastore_dir}, should "
            + ["not", "still"][retains_datastore]
            + "after the removal operation"
        )


@pytest.mark.parametrize(
    "rm_approach", ["rm-implicit-vdata", "rm-explicit-vdata", "rm-data-store"]
)
def test_fetch_and_remove(dummy_file_repo, rm_approach, cli_app):
    # in this test, the fetch operation is used to fetch all files in the registry
    # and we vary the rm approach

    version = "1.0"
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary", version)
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")
    with custom_datadir(data_dir):
        # fetch the data
        app.fetch(dummy_file_repo.src_file_dir.primary, expect_success=True)
        check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        # confirm that if we call fetch again, that we still consider it a success
        app.fetch(dummy_file_repo.src_file_dir.primary, expect_success=True)
        check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        # now, lets remove the version-data
        # -> currently, we explicitly test a handful of scenarios where we expect things
        #    to fail. It might be nice to split these cases out into separate tests
        if "rm-data-store" == rm_approach:
            rm_method, rm_args, retains_datastore = app.rm_datastore, (), False
        elif "rm-implicit-vdata" == rm_approach:
            rm_method, rm_args, retains_datastore = app.rm_vdata, (None,), True
        elif "rm-explicit-vdata" == rm_approach:
            rm_method, rm_args, retains_datastore = app.rm_vdata, (version,), True

            # extra scenario worth testing - trying to remove values that don't exist
            app.rm_vdata(version="9999999999999.0", expect_success=False)
            check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)
        else:
            raise RuntimeError("unexpected rm_approach")

        # we expect the remove command to report success, when we omit the force flag,
        # but not actually do anything
        rm_method(*rm_args, omit_force=True, expect_success=True)
        check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        # we expect the following case to fail
        with open(_get_lockfile_path(data_dir), "w") as f:
            f.write("a dummy lockfile")
        rm_method(*rm_args, expect_success=False)
        check_version_data_dir(
            data_dir, version, lockfile_should_exist=True, file_set=_DUMMY_SET_PRIMARY
        )
        os.remove(_get_lockfile_path(data_dir))

        # now we expect it to succeed
        rm_method(*rm_args, expect_success=True)

        _check_removal(data_dir, version, retains_datastore)


@pytest.mark.parametrize(
    "fetch_subset",
    [
        "fetch-single",
        "fetch-single-then-all",
        "fetch-all-then-single",
        "fetch-all-explicit-list",
    ],
)
def test_fetch_subset_and_remove(dummy_file_repo, fetch_subset, cli_app):
    # in this test, the fetch operation is used to fetch a named subset of files
    # in the registry and we use a single rm approach
    version = "1.0"
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary", version)
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")
    with custom_datadir(data_dir):
        # fetch the data
        if fetch_subset == "fetch-single":
            # call it twice to ensure we consider the operation a success each time
            for i in range(2):
                app.fetch(
                    dummy_file_repo.src_file_dir.primary,
                    file_list=[_DUMMY_SET_PRIMARY[1][0]],
                    expect_success=True,
                )
                check_version_data_dir(
                    data_dir, version, file_set=[_DUMMY_SET_PRIMARY[1]]
                )

        elif fetch_subset == "fetch-single-then-all":
            app.fetch(
                dummy_file_repo.src_file_dir.primary,
                file_list=[_DUMMY_SET_PRIMARY[1][0]],
                expect_success=True,
            )
            check_version_data_dir(data_dir, version, file_set=[_DUMMY_SET_PRIMARY[1]])

            # fetch all data
            app.fetch(dummy_file_repo.src_file_dir.primary, expect_success=True)
            check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        elif fetch_subset == "fetch-all-then-single":
            # fetch all data
            app.fetch(dummy_file_repo.src_file_dir.primary, expect_success=True)
            check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

            # fetch a single data file
            app.fetch(
                dummy_file_repo.src_file_dir.primary,
                file_list=[_DUMMY_SET_PRIMARY[1][0]],
                expect_success=True,
            )
            check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        elif fetch_subset == "fetch-all-explicit-list":
            # fetch all data files
            app.fetch(
                dummy_file_repo.src_file_dir.primary,
                file_list=[pair[0] for pair in _DUMMY_SET_PRIMARY],
                expect_success=True,
            )
            check_version_data_dir(data_dir, version, file_set=_DUMMY_SET_PRIMARY)

        else:
            raise RuntimeError("unexpected fetch_subset")

        app.rm_vdata(None, expect_success=True)
        _check_removal(data_dir, version, retains_datastore=True)


@pytest.mark.parametrize(
    "src_file_dir_key",
    [
        pytest.param("replace", id="cksum-mismatch"),
        pytest.param("rename", id="no-src-file"),
    ],
)
def test_fetch_fail(src_file_dir_key, dummy_file_repo, cli_app):
    # here we intentionally use a file registry and a mismatched source-directory (where
    # files are fetched from). Essentially we want to ensure correct (and graceful)
    # behavior in 2 failure modes:
    # 1. somehow the checksum is mismatched. This might happen if a file got corrupted in
    #    a download OR (more likely) we made an error while creating the registry file.
    # 2. somwhow the file can't be fetched. This might happen for a range of reasons
    #    such as internet connectivity issues or server issues. Alternatively it could
    #    happen if we make a mistake (e.g. update the registry, but forget to upload the
    #    data file)
    nominal_version = "1.0"
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary", nominal_version)
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")
    with custom_datadir(data_dir):
        success = app.fetch(dummy_file_repo.src_file_dir.get(src_file_dir_key))
    assert not success

    # the directories should exist, but the file should not be "downloaded". We also
    # confirm that there isn't a lockfile (i.e. we exit gracefully)
    check_version_data_dir(data_dir, nominal_version, lockfile_should_exist=False)
    assert not os.path.isfile(
        _get_managed_file(data_dir, nominal_version, _DUMMY_SET_REPLACE[1][0])
    )


def test_fetch_fail_locked(dummy_file_repo, cli_app):
    # confirm that the fetch command will fail if a lockfile exists
    nominal_version = "1.0"
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary", nominal_version)
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")
    # let's create a lockfile
    os.makedirs(data_dir)
    with open(_get_lockfile_path(data_dir), "w") as f:
        f.write("a dummy lockfile")

    with custom_datadir(data_dir):
        success = app.fetch(dummy_file_repo.src_file_dir.primary)
    assert not success, "Failure is expected when a lockfile exists"
    if os.path.isdir(_get_version_dir(data_dir, nominal_version)):
        raise AssertionError(
            "the tool should not create a version directory when a lock file exists"
        )


def test_fetch_untracked(dummy_file_repo, cli_app):
    # test the fetch operation is used to fetch files to an untracked directory
    version = "1.0"
    app = dummy_file_repo.cli_app_with_overrides(cli_app, "primary", version)
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")
    with custom_datadir(data_dir):
        # first we download everything
        dest_dir_1 = os.path.join(dummy_file_repo.test_dir, "my-untracked-dir-all")
        for i in range(2):
            app.fetch(
                dummy_file_repo.src_file_dir.primary,
                untracked_dest_dir=dest_dir_1,
                expect_success=True,
            )
            assert not os.path.isfile(data_dir)
            check_version_data_dir_contents(dest_dir_1, file_set=_DUMMY_SET_PRIMARY)

        # now confirm that we can download a subset
        dest_dir_2 = os.path.join(dummy_file_repo.test_dir, "my-untracked-dir-single")
        for i in range(2):
            app.fetch(
                dummy_file_repo.src_file_dir.primary,
                untracked_dest_dir=dest_dir_2,
                file_list=[_DUMMY_SET_PRIMARY[1][0]],
                expect_success=True,
            )
            assert not os.path.isfile(data_dir)
            check_version_data_dir_contents(
                dest_dir_2, file_set=[_DUMMY_SET_PRIMARY[1]]
            )


def is_linked(*paths):
    if len(paths) < 2:
        raise TypeError("is_linked() must have at least 2 arguments")
    try:
        # if ref has a value of 0, then we can't actually make a meaningful
        # comparison. There isn't any obvious behavior in this scenario
        ref = os.stat(paths[0], follow_symlinks=True).st_ino
        return all(
            ref == os.stat(path, follow_symlinks=True).st_ino for path in paths[1:]
        )
    except FileNotFoundError:
        return False


def test_multiversion(dummy_file_repo, cli_app):
    # test what happens when we fetch multiple sets of files
    # - in the future, it might be nice to break this test up into smaller pieces

    # from a realism perspective, primary -> rename -> replace may make more sense, but
    # the current order seems more likely to catch an error
    version_kind_map = {"1.0": "primary", "2.0": "replace", "3.0": "rename"}
    versions = tuple(version_kind_map.keys())
    app_v1, app_v2, app_v3 = [
        dummy_file_repo.cli_app_with_overrides(cli_app, kind, version)
        for version, kind in version_kind_map.items()
    ]
    data_dir = os.path.join(dummy_file_repo.test_dir, "my-data-dir")

    def _basic_datastore_check(expected_versions, last_op, last_op_version):
        err_msg = (
            f"This check is performed right after performing the `{last_op}` "
            "operation with grdata, specialized for simulated version "
            f"{last_op_version} of grackle. (That version is associated with the "
            f"{version_kind_map[last_op_version]!r} dummy fileset)"
        )

        for ver in expected_versions:
            check_version_data_dir(
                data_dir,
                ver,
                file_set=_DUMMY_SET_TUPLE.get(version_kind_map[ver]),
                err_msg=err_msg,
            )
        for ver in version_kind_map.keys():
            if ver not in expected_versions:
                ver_dir = _get_version_dir(data_dir, ver)
                if os.path.isdir(ver_dir):
                    raise AssertionError(
                        f"{err_msg}\n\n"
                        "The version-directorty, {ver_dir}, should not exist!"
                    )

    with custom_datadir(data_dir):
        # step 1: load data associated with v1 (the `primary` fileset)
        app_v1.fetch(dummy_file_repo.src_file_dir.primary, expect_success=True)
        _basic_datastore_check(versions[:1], "fetch", last_op_version=versions[0])

        # step 2: load data associated with v2 (the `replace` fileset)
        app_v2.fetch(dummy_file_repo.src_file_dir.replace, expect_success=True)
        _basic_datastore_check(versions[:2], "fetch", last_op_version=versions[1])

        # confirming linking... (it might be better to check disk-usage and be more
        # agnostic about deduplication)
        assert is_linked(
            _get_managed_file(data_dir, versions[0], _DUMMY_SET_TUPLE.primary[0][0]),
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[0][0]),
        ), "the file-0.txt files should all be linked"
        assert not is_linked(
            _get_managed_file(data_dir, versions[0], _DUMMY_SET_TUPLE.primary[1][0]),
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[1][0]),
        ), "the file-1.txt files are expected to hold different contents"

        # step 3: load data associated with v3 (the `rename` fileset)
        app_v3.fetch(dummy_file_repo.src_file_dir.rename, expect_success=True)
        _basic_datastore_check(versions, "fetch", last_op_version=versions[2])
        # checking linking... (it might be better to check disk-usage and be more
        # agnostic about the fact that we use linking for deduplication)
        assert is_linked(
            _get_managed_file(data_dir, versions[0], _DUMMY_SET_TUPLE.primary[0][0]),
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[0][0]),
            _get_managed_file(data_dir, versions[2], _DUMMY_SET_TUPLE.rename[0][0]),
        ), "the file-0.txt files should all be linked"
        assert is_linked(
            _get_managed_file(data_dir, versions[0], _DUMMY_SET_TUPLE.primary[1][0]),
            _get_managed_file(data_dir, versions[2], _DUMMY_SET_TUPLE.rename[1][0]),
        ), "the file-1.txt files from primary and rename filesets should be linked"
        assert not is_linked(
            _get_managed_file(data_dir, versions[0], _DUMMY_SET_TUPLE.primary[1][0]),
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[1][0]),
        ), "file-1.txt file from the replace fileset should not be linked to anything"

        # step 4: remove data associated with v1 (the `primary` fileset)
        # -> we EXPLICITLY use a different app version to remove this data
        app_v3.rm_vdata(versions[0], expect_success=True)
        _basic_datastore_check(
            [versions[1], versions[2]], "rm-vdata", last_op_version=versions[0]
        )
        assert is_linked(
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[0][0]),
            _get_managed_file(data_dir, versions[2], _DUMMY_SET_TUPLE.rename[0][0]),
        ), "remaining file-0.txt files should remain linked"
        assert not is_linked(
            _get_managed_file(data_dir, versions[1], _DUMMY_SET_TUPLE.replace[1][0]),
            _get_managed_file(data_dir, versions[2], _DUMMY_SET_TUPLE.rename[1][0]),
        ), "remaining file-1.txt files should remain unlinked"

        # step 5: remove data associated with v3 (the `rename` fileset)
        # -> we EXPLICITLY use a different app version to remove this data
        app_v2.rm_vdata(versions[2], expect_success=True)
        _basic_datastore_check([versions[1]], "rm-vdata", last_op_version=versions[2])
