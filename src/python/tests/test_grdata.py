########################################################################
#
# Test the grdata command line tool for fetching test files
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
import dataclasses
import filecmp
import functools
import hashlib
import io
import operator
import os
import shutil
import subprocess
import sys
from textwrap import indent
from typing import Callable, Collection, ClassVar, Dict, List, Optional, Tuple

import pytest

# a goal here is to be able to run this test without installing pygrackle!
# -> in the near future, we will install grdata as a standalone command-line
#    script & it might be nice to test script without installing pygrackle
# -> when that time comes, we will modify the logic within the nominal_app fixture
try:
    from pygrackle.__main__ import main as grdata_main
    from pygrackle.utilities.grdata import _parse_file_registry
except ImportError:
    grdata_main, _parse_file_registry = None, None


# _ENV_VAR holds environment variables that can influence data directory location
_ENV_VARS = ("GRACKLE_DATA_DIR", "XDG_DATA_HOME", "HOME")
_CKSUM_ALG = "sha1"


@contextlib.contextmanager
def modified_env(new_env_vals, extra_cleared_variables=None):
    """
    Temporarily overwrite the environment variables. This is necessary to test C
    extensions that rely upon the environment variables
    """

    _UNSET = object()
    modify_set = set(new_env_vals.keys())
    if extra_cleared_variables is not None:
        modify_set.update(extra_cleared_variables)

    original_vals = {var: os.environ.pop(var, _UNSET) for var in modify_set}
    for var, value in new_env_vals.items():
        os.environ[var] = value
    try:
        yield
    finally:
        for var, value in original_vals.items():
            if value is _UNSET:
                del os.environ[var]
            else:
                os.environ[var] = value


@dataclasses.dataclass(frozen=True)
class CLIApp:
    """A wrapper around the grdata tool"""

    cli_launcher_args: Tuple[str, ...]
    # the following are all used to override the choices within grdata
    override_datadir: Optional[str] = None
    override_grackle_version: Optional[str] = None
    override_file_registry_file: Optional[str] = None
    override_base_url: Optional[str] = None

    # maps flags to override-attrs
    _OVERRIDE_FLAG_MAP: ClassVar[Dict[str, str]] = {
        "--override-base-url": "override_base_url",
        "--testing-override-file-registry-file": "override_file_registry_file",
        "--testing-override-grackle-version": "override_grackle_version",
    }

    def all_overriden(self) -> bool:
        attrs = list(self._OVERRIDE_FLAG_MAP.values())
        attrs.append("override_datadir")
        return all(getattr(self, attr) is not None for attr in attrs)

    def replace(self, **kwargs) -> "CLIApp":
        """return a new CLIApp instance replacing specified fields with new values"""
        return dataclasses.replace(self, **kwargs)

    @functools.cached_property
    def _common_args(self) -> List[str]:
        args = [elem for elem in self.cli_launcher_args]
        for flag, attr in self._OVERRIDE_FLAG_MAP.items():
            val = getattr(self, attr)
            if isinstance(val, str):
                args.extend([flag, val])
            elif val is not None:  # this is a shockingly common mistake
                raise TypeError(f"the {attr} attribute must be a string")
        return args

    def __call__(self, *subcommand_args: str, expect_success: bool = True) -> str:
        """Invoke a subcommand"""
        __tracebackhide__ = True  # have pytest hide noisy tracebacks

        if len(subcommand_args) == 0:
            name = self.__class__.__name__
            raise ValueError(f"{name}.__call__(...) requires at least 1 positional arg")
        else:
            all_args = self._common_args.copy()
            all_args.extend(subcommand_args)

        with contextlib.ExitStack() as stack:
            if self.override_datadir is not None:
                new_env_vals = {"GRACKLE_DATA_DIR": self.override_datadir}
                stack.enter_context(modified_env(new_env_vals))

            if len(self.cli_launcher_args) > 0:
                tmp = subprocess.run(all_args, capture_output=True)
                returncode = tmp.returncode
                stdout, stderr = tmp.stdout.decode("ascii"), tmp.stderr.decode("ascii")
            elif grdata_main is None:
                raise RuntimeError("can't be configured to run without subprocess")
            else:
                f_out, f_err = [
                    stack.enter_context(contextlib.redirect_stdout(io.StringIO())),
                    stack.enter_context(contextlib.redirect_stderr(io.StringIO())),
                ]
                returncode = grdata_main(all_args)
                stdout, stderr = f_out.getvalue(), f_err.getvalue()

            if (returncode == 0) == expect_success:
                return stdout.rstrip()

            # handle error message
            details = {"args": all_args}
            for v in _ENV_VARS:
                details[f"environ[{v!r}]"] = os.environ.get(v, "<unset>")
            details.update({"code": returncode, "stdout": stdout, "stderr": stderr})

            outcome = "succeeded" if returncode == 0 else "failed"
            msg_lines = [f"grdata unexpectedly {outcome}:\n"]
            for key, val in details.items():
                multiline_str = isinstance(val, str) and "\n" in val
                s = f"\n{indent(val, '   > ')}" if multiline_str else repr(val)
                msg_lines.append(f"  {key}: {s}")
            raise RuntimeError("\n".join(msg_lines))


@dataclasses.dataclass
class DummyFileSpec:
    contents_str: str
    cksum: str


@functools.lru_cache
def _spec_map(return_v1: bool, cksum_alg: str) -> Dict[str, DummyFileSpec]:
    def make_spec(variant: str, trailing_text: str = ""):
        content = f"I am a test-file: Variant {variant}\n{trailing_text}"
        content_bytes = content.encode("ascii")

        if shutil.which("shasum") is not None:
            _alg = {"sha1": "1", "sha256": "256"}[cksum_alg]
            args = ["shasum", "--algorithm", _alg, "-"]
            kw = {"input": content_bytes, "check": True, "capture_output": True}
            rslt_str = subprocess.run(args, **kw).stdout.rstrip().decode("utf8")
            if not rslt_str.endswith("  -"):
                raise RuntimeError(f"output of shasum was unexpected: '{rslt_str}'")
            cksum = f"{cksum_alg}:{rslt_str[:-3].lower()}"
        else:
            hash_obj = hashlib.new(cksum_alg)
            hash_obj.update(content_bytes)
            cksum = f"{cksum_alg}:{hash_obj.hexdigest()}"
        return DummyFileSpec(content, cksum)

    fA_spec = make_spec("A")
    fB_spec = make_spec("B")
    if return_v1:
        return {"A.txt": fA_spec, "B.txt": fB_spec}
    else:
        mut_fB_spec = make_spec(1, trailing_text="version 2 of file\n")
        return {"A.txt": fA_spec, "old-B.txt": fB_spec, "B.txt": mut_fB_spec}


def _setup_cliapp(root_path: os.PathLike, return_v1: bool, prototype: CLIApp):
    """sets up an CLIApp instance to work with dummy files

    This is intended to be used by a fixture. In more detail, this function
    first sets up some directories within root_path:
      - the fixture_dir isn't supposed to be mutated by tests. It contains a
        registry-file & a directory containing files from a dummy-file-set.
        (the choice of dummy file-set is controlled by the `return_v1` arg)
      - the test-dir is reserved for mutation by the tests.
    After finishing this, the function returns an instance of CLIApp that is
    properly configured to use this information. This instance is copied from
    the instance passed as the `prototype` argument.

    This function is callable twice with the same `root_path` argument as long
    as `return_v1` has different values each time. It is instructive to look at
    the following sketch of the directory structure resulting from this scenario:

     root_path/
      ├── fixture_dir/                # <- this directory holds inputs for tests
      │    ├── v1_file_registry.txt      # registry for "v1" dummy-file-set
      │    ├── v1-ref-files/             # holds files in "v1" dummy-file-set
      │    │    ├── A.txt
      │    │    └── B.txt
      │    ├── v2_file_registry.txt      # registry for "v2" dummy-file-set
      │    └── v2-ref-files/             # holds files in "v2" dummy-file-set
      │         ├── A.txt                # matches contents of v1-ref-files/A.txt
      │         ├── old-B.txt            # matches contents of v1-ref-files/B.txt
      │         └── B.txt                # unique contents!
      └── test-dir/                   # <- tests can mutate this directory
    """

    test_dir = os.path.abspath(os.path.join(root_path, "test-dir"))
    data_dir = os.path.abspath(os.path.join(test_dir, "my-data-dir"))
    fixture_dir = os.path.abspath(os.path.join(root_path, "fixture_dir"))

    if not os.path.isdir(test_dir):
        os.mkdir(test_dir)
        os.mkdir(fixture_dir)

    spec_map = _spec_map(return_v1=return_v1, cksum_alg=_CKSUM_ALG)
    vstr = "v1.0" if return_v1 else "v2.0"
    src_file_dir = os.path.join(fixture_dir, f"{vstr}-ref-files")
    registry_path = os.path.join(fixture_dir, f"{vstr}_file_registry.txt")

    os.mkdir(src_file_dir)
    registry_lines = []
    for fname, file_spec in spec_map.items():
        with open(os.path.join(src_file_dir, fname), "w") as f:
            f.write(file_spec.contents_str)
        registry_lines.append(f'{{"{fname}", "{file_spec.cksum}"}}')
    with open(registry_path, "w") as f:
        print(*registry_lines, sep=",\n", end="\n", file=f)

    return prototype.replace(
        override_datadir=os.fsdecode(data_dir),
        override_grackle_version=vstr,
        override_file_registry_file=os.fsdecode(registry_path),
        override_base_url=f"file:/{os.fsdecode(src_file_dir)}",
    )


@pytest.fixture(scope="module")
def nominal_app():
    if True:
        assert sys.executable is not None
        return CLIApp(cli_launcher_args=(sys.executable, "-m", "pygrackle"))
    else:  # this branch makes the tests run ~100 times faster but is less general
        return CLIApp(cli_launcher_args=())


@pytest.fixture
def app_v1(tmp_path, nominal_app):
    cli_app = _setup_cliapp(root_path=tmp_path, return_v1=True, prototype=nominal_app)
    yield cli_app


@pytest.fixture
def app_v2(tmp_path, nominal_app):
    cli_app = _setup_cliapp(root_path=tmp_path, return_v1=False, prototype=nominal_app)
    yield cli_app


def test_checksum_kind(nominal_app):
    # we are confirming that the hardcoded checksum in grdata hasn't changed
    cksum_kind = nominal_app("show-checksum-kind")
    assert cksum_kind == _CKSUM_ALG


def test_overrides(app_v1):
    # essentially, we are checking that the overrides all work
    if not app_v1.all_overriden():
        pytest.skip()

    grackle_version_str = app_v1("show-grackle-version")
    assert grackle_version_str.strip() == app_v1.override_grackle_version

    datadir = app_v1("show-data-dir")
    assert os.path.abspath(app_v1.override_datadir) == os.path.abspath(datadir)

    full_registry_str = app_v1("show-known-reg")
    with open(app_v1.override_file_registry_file, "r") as f:
        ref_full_registry_str = f.read().rstrip()
    assert full_registry_str.rstrip() == ref_full_registry_str


def test_calcreg(app_v1):
    # we are confirming that calling the calc-reg subcommand will produce a registry
    # consistent with the known hashes

    assert app_v1.override_base_url is not None
    if not app_v1.override_base_url.startswith("file:"):
        pytest.skip()
    repository_path = app_v1.override_base_url[5:]

    registry_str = app_v1("calc-reg", f"--hash-alg={_CKSUM_ALG}", repository_path)
    registry = _parse_file_registry(io.StringIO(registry_str))
    ref_registry = _parse_file_registry(app_v1.override_file_registry_file)

    expected_fname_set = set(ref_registry.keys())
    if len(expected_fname_set.symmetric_difference(registry.keys())) != 0:
        raise AssertionError(
            "calc-reg subcommand produced a registry with entries for the dummy files: "
            f"{sorted(registry)}. Expected the registry to have entries for "
            f"{sorted(expected_fname_set)}."
        )

    for fname, cksum in registry.items():
        ref = ref_registry[fname]
        if cksum != ref:
            raise AssertionError(
                f"grdata computed an unexpected checksum for the dummy-file, {fname}\n"
                f"  expected: {ref}\n  actual: {cksum}"
            )


def _get_version_dir(app: CLIApp) -> str:
    args = (app.override_datadir, "data-store-v1", app.override_grackle_version)
    return os.path.join(*args)


def check_fetched_dir_contents(
    app: CLIApp,
    fetch_dir: Optional[str] = None,
    *,
    fname_subset: Optional[Collection[str]] = None,
    errmsg_writer: Optional[Callable] = None,
):
    __tracebackhide__ = True  # suppress noisy pytest tracebacks

    errmsg_writer = AssertionError if errmsg_writer is None else errmsg_writer
    if not (app.all_overriden() and app.override_base_url.startswith("file:")):
        raise RuntimeError("app argument doesn't satisfy preconditions")
    ref_dir = app.override_base_url[5:]
    fetch_dir = _get_version_dir(app) if fetch_dir is None else fetch_dir
    fnames = os.listdir(ref_dir) if fname_subset is None else fname_subset

    for fname in fnames:
        path, ref_path = os.path.join(fetch_dir, fname), os.path.join(ref_dir, fname)
        if not os.path.isfile(ref_path):
            raise RuntimeError(f"reference file is missing: {ref_path}")
        elif not os.path.isfile(path):
            raise errmsg_writer(f"{path}, doesn't exist after invoking grdata")
        elif not filecmp.cmp(path, ref_path, shallow=False):
            raise errmsg_writer(f"{path}, doesn't have the correct contents")

    if (fname_subset is None) and (len(fnames) != len(os.listdir(fetch_dir))):
        raise errmsg_writer(f"dir doesn't contain correct entries: {fetch_dir}")


def check_version_data_dir(
    app: CLIApp,
    fname_subset: Optional[Collection[str]] = None,
    err_msg: Optional[str] = None,
):
    __tracebackhide__ = True  # suppress noisy pytest tracebacks
    if not app.all_overriden():
        raise RuntimeError("app argument doesn't satisfy preconditions")

    required_paths = [
        ("data directory", app.override_datadir),
        ("user-data directory", os.path.join(app.override_datadir, "user-data")),
        ("version data directory", _get_version_dir(app)),
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

    check_fetched_dir_contents(
        app, fname_subset=fname_subset, errmsg_writer=prep_assertion_err
    )


def test_fetch(app_v1):
    app_v1("fetch")  # fetch all files (v1.0/A.txt & v1.0/B.txt)
    check_version_data_dir(app_v1)
    app_v1("fetch")  # does nothing (but should still succeed)
    app_v1("fetch", "A.txt")  # does nothing (but should still succeed)
    check_version_data_dir(app_v1)  # <- sanity check!


def test_fetch_single(app_v1):
    app_v1("fetch", "A.txt")  # fetch v1.0/A.txt
    check_version_data_dir(app_v1, fname_subset=["A.txt"])
    app_v1("fetch", "B.txt")  # fetch v1.0/B.txt
    check_version_data_dir(app_v1)

    app_v1("fetch", "A.txt")  # does nothing (but should still succeed)
    app_v1("fetch")  # does nothing (but should still succeed)
    check_version_data_dir(app_v1)  # <- sanity check!


def test_fetch_all_explicit(app_v1):
    # we explicitly list all of the files
    app_v1("fetch", "A.txt", "B.txt")
    check_version_data_dir(app_v1)
    app_v1("fetch")  # does nothing (but should still succeed)
    check_version_data_dir(app_v1)  # <- sanity check!


def test_fetch_fail_unknown_file(app_v1):
    app_v1("fetch", "not-a-file", expect_success=False)  # <- should fail


def test_fetch_fail_cksum_mismatch(app_v1, app_v2):
    # create a new application where file-registry and base-url are mismatched
    app_mismatch = app_v1.replace(
        override_grackle_version="mismatched-registry",
        override_file_registry_file=app_v2.override_file_registry_file,
    )

    # the following will work since A.txt is identical for app_v1 & app_v2
    app_mismatch("fetch", "A.txt")
    check_version_data_dir(app_mismatch, fname_subset=["A.txt"])

    # the following won't work (B.txt has different checksums for app_v1 & app_v2)
    app_mismatch("fetch", "B.txt", expect_success=False)
    check_version_data_dir(app_mismatch, fname_subset=["A.txt"])


def test_fetch_untracked_all(tmp_path, app_v1):
    dest = os.fsdecode(tmp_path / "untracked-dir-all")
    assert not os.path.isdir(dest)
    for _ in range(2):  # the second time should do nothing (but should succeed)
        app_v1("fetch", "--untracked-dest-dir", dest)
        check_fetched_dir_contents(app_v1, dest)


def test_fetch_untracked_single(tmp_path, app_v1):
    dest = os.fsdecode(tmp_path / "untracked-dir-single")
    assert not os.path.isdir(dest)
    for i in range(2):  # the second time should do nothing (but should succeed)
        app_v1("fetch", "--untracked-dest-dir", dest, "A.txt")
        check_fetched_dir_contents(app_v1, dest, fname_subset=["A.txt"])


def test_deduplication(app_v1, app_v2):
    basedir = os.path.join(app_v1.override_datadir, "data-store-v1")

    def is_linked(path_a, path_b):
        path_a, path_b = os.path.join(basedir, path_a), os.path.join(basedir, path_b)
        ref = os.stat(path_a, follow_symlinks=True).st_ino
        return ref == os.stat(path_b, follow_symlinks=True).st_ino

    # setup app_v10, app_v11, & app_v20 (all have same registry but app_v20)
    app_v10, app_v20 = app_v1, app_v2
    app_v11 = app_v1.replace(override_grackle_version="v1.1")

    app_v10("fetch", "B.txt")  # setup v1.0/A.txt
    check_version_data_dir(app_v10, ["B.txt"])
    app_v11("fetch")  # setup v1.1/A.txt & v1.1/B.txt
    check_version_data_dir(app_v11)
    assert is_linked("v1.0/B.txt", "v1.1/B.txt")

    app_v20("fetch")  # setup v2.0/A.txt, v2.0/old-B.txt, & v2.0/B.txt
    check_version_data_dir(app_v20)
    assert is_linked("v1.1/A.txt", "v2.0/A.txt")
    assert not is_linked("v1.1/B.txt", "v2.0/B.txt")
    assert not is_linked("v1.1/B.txt", "v2.0/old-B.txt")

    app_v10("fetch", "A.txt")  # setup v1.0/B.txt
    check_version_data_dir(app_v10)
    assert is_linked("v1.0/A.txt", "v1.1/A.txt")
    assert is_linked("v1.0/A.txt", "v2.0/A.txt")
