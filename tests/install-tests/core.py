"""
Define the core testing logic
"""

# to maximize portability, only use built-in modules present in python 3.6.1 and newer
from collections import ChainMap
from enum import Enum
import fnmatch
import functools
import os
import re
import shlex
import sys
from typing import (
    Callable,
    Dict,
    Iterator,
    Mapping,
    NamedTuple,
    Optional,
    Tuple,
    Union,
)

from docker import DockerContainer, GRACKLE_ROOT, import_standalone_module

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ImportError:
        import vendored_tomli as tomllib


common = import_standalone_module(
    dir_path=os.path.join(GRACKLE_ROOT, "scripts", "ci"), module_name="common"
)
ScriptError = common.ScriptError
logger = common.logger
make_cmd_summary_str = common.make_cmd_summary_str


# constants that need to remain synchronized with Dockerfile
#
# -> in an ideal world, we would move away from hardcoding all these things
#    -> we could infer DockerImage members by scanning Dockerfile
#    -> we could create a json file as part of the Dockerfile that holds things like:
#       -> _LOCAL_LIB
#       -> Grackle_DIR (not defined here, but its hardcoded in the conf files)

_USERNAME = "gr-user"
_LOCAL_LIB = f"/home/{_USERNAME}/local/lib"


class DockerImage(Enum):
    """Constants describing the docker images"""

    baseline = "baseline"
    shared_build = "shared_build"
    shared_install = "shared_install"
    shared_static_install = "shared_static_install"
    static_build = "static_build"
    static_install = "static_install"
    static_shared_install = "static_shared_install"
    static_install_conftimeprefix = "static_install_conftimeprefix"

    def layer_name(self):
        # the layer name in ./Dockerfile corresponding to this image
        return self._value_

    def tag_name(self):
        return f"grackle_install_test_{self._value_}"

    def __repr__(self):
        return f"<{self.__class__.__name__}.{self._value_}>"

    def get_grackle_dir(self):
        return f"/home/{_USERNAME}/grackle"


class TestItem(NamedTuple):
    """Represents a test item"""

    image: DockerImage

    # to construct the sample project, we copy all paths in
    # [f"{GRACKLE_ROOT}/{e} for in src_files]
    src_files: Tuple[str, ...]

    # the collection of commands to execute for the build
    cmds: Tuple[str, ...]

    # specifies the path to the output binary relative to the root of the source dir.
    # A value of None, indicates that a build is expected to fail
    output_binary: Optional[str]

    # specifies whether LD_LIBRARY_PATH needs to be overwritten
    override_LD_LIBRARY_PATH: bool = False

    @property
    def expect_success(self):
        return self.output_binary is not None


# we are going to be fairly pedantic about toml parsing
# define 2 functions to check the data types specified in the toml file
def _is_table_arr(v) -> bool:
    # does v correspond to a TOML array of tables?
    # -> TOML doesn't support an empty array of tables
    # -> if 1 array entry is a table then all array entries are tables
    return isinstance(v, list) and len(v) > 0 and isinstance(v[0], dict)


def _is_str_arr(v) -> bool:
    return isinstance(v, list) and all(isinstance(elem, str) for elem in v)


# lists the requirements for all allowed keys in one of the TOML test tables
# (and optionally specify a coercion function)
_REQS = {
    "image": ((str, DockerImage), DockerImage, "a string"),
    "cmds": (_is_str_arr, tuple, "an array of strings"),
    "src_file": (str, None, "a string"),
    "output_binary": (str, None, "a string"),
    "override_LD_LIBRARY_PATH": (bool, None, "a boolean"),
}


def _coerce_toml_table(table: Mapping, err_prefix=""):
    out = {}
    for key, value in table.items():
        try:
            classinfo, coerce_fn, descr = _REQS[key]
        except KeyError:
            raise ValueError(f"{err_prefix}unexpected key: {key}") from None

        if isinstance(classinfo, (tuple, type)):
            has_right_type = isinstance(value, classinfo)
        else:
            has_right_type = classinfo(value)
        if not has_right_type:
            raise ValueError(f"{err_prefix}{key} must be {descr}")
        out[key] = value if coerce_fn is None else coerce_fn(value)
    return out


def test_cases_from_conf(path: str) -> Iterator[Tuple[str, TestItem]]:
    """Constructs TestItem from the specified configuration toml file"""

    # before we actually parse the config file, let's build a list of all other
    # files in the the config file
    # -> we assume that these are parts of the sample source directory
    # -> we store the paths to these files so that they're relative to GRACKLE_ROOT
    abs_conf_dir_path = os.path.dirname(os.path.abspath(path))
    _rel_conf_dir_path = os.path.relpath(abs_conf_dir_path, start=GRACKLE_ROOT)
    common_src_files = []
    with os.scandir(abs_conf_dir_path) as it:
        for entry in it:
            if entry.is_dir(follow_symlinks=True):
                raise ValueError(f"dir holding {path} can't also hold subdirectories")
            elif os.path.samefile(entry.path, path):
                continue  # <- we won't need to copy the conf file itself
            else:
                common_src_files.append(os.path.join(_rel_conf_dir_path, entry.name))

    with open(path, "r") as f:
        raw = tomllib.load(f)

    err_prefix = f"Toml File Error {path}: "

    # read top-level entities and do error-checks
    common_table = raw.get("common", {})
    if not isinstance(common_table, dict):
        raise ValueError(f'{err_prefix}"common" key must provide a table, if specified')
    common_table = _coerce_toml_table(common_table, err_prefix=err_prefix)

    try:
        test_arr = raw["test"]
    except KeyError:
        raise ValueError(f'{err_prefix} "test" key is missing') from None
    if not _is_table_arr(test_arr):
        raise ValueError(
            f'{err_prefix}"test" key must correspond to an array of tables'
        )

    name_prefix = os.path.basename(abs_conf_dir_path)
    for table in test_arr:
        table = _coerce_toml_table(table, err_prefix=err_prefix)
        full_table = ChainMap(table, common_table, TestItem._field_defaults)
        if "src_file" not in full_table:
            raise ValueError(f'{err_prefix}"src_file" wasn\'t specified')

        kwargs = {k: v for k, v in full_table.items() if k != "src_file"}
        kwargs["src_files"] = tuple(
            [os.path.join("src", "example", full_table["src_file"])] + common_src_files
        )
        test_item = TestItem(**kwargs)

        test_name = f"{name_prefix}.{test_item.image._value_}"
        yield (test_name, test_item)


class TestResult(NamedTuple):
    action_during_failure: Optional[str]
    issue: Union[BaseException, str, None]

    @classmethod
    def Success(cls) -> "TestResult":
        return cls(None, None)

    @classmethod
    def Fail(
        cls, action_during_failure: str, issue: Union[BaseException, str, None] = None
    ) -> "TestResult":
        return cls(action_during_failure, issue)

    @classmethod
    def ProgrammingError(cls, exception: BaseException) -> "TestResult":
        return cls(None, exception)

    def is_success(self):
        return self.action_during_failure is None and self.issue is None


def _is_quoted(v: str, quote: str):
    return v[:1] == v[:-1] == quote


def _exec_container_cmd(
    *args,
    container: DockerContainer,
    print_info: Callable = print,
    expect_success: bool = True,
    **call_kwargs,
) -> TestResult:
    """
    Helper function to run_test that appropriately dispatches to
    DockerContainer.call
    """
    # ugh, we have like 3 different functions that do this kind of thing...
    # maybe we should refactor?
    #
    # Currently, this is necessary to make sure we package up a TestResult

    if "log" in call_kwargs or "capture_output" in call_kwargs:
        raise RuntimeError("log and capture_output shouldn't be in call_kwargs")
    _msg = make_cmd_summary_str(
        args=args,
        include_outer_env=call_kwargs.get("include_outer_env", True),
        env=call_kwargs.get("env", None),
        cwd=call_kwargs.get("cwd", container.home_dir),
    )

    # (I'm not sure we really want to support this...)
    # we support bash-style inline env var assignment, where we write something like
    #   `VAR=VAL <cmd> <args...>` or `VAR1=VAL1 VAR2=VAL2 VAR3=VAL3 <cmd> <args...>`
    # Importantly, bash
    #   - does all variable substitution before evaluating a line
    #   - takes the last definition if a Variable name is defined multiple times
    overrides = {}
    while len(args) > 1:
        match = re.match(r"^([a-zA-Z0-9_]+)=(.*)$", args[0])
        if match is None:
            break
        _var, _val = match.groups()
        if _is_quoted(_val, "'"):
            overrides[_var] = _val[2:-1]
        elif "$" not in _val and "\\" not in _val:
            overrides[_var] = _val[2:-1] if _is_quoted(_val, '"') else _val
        else:
            _echo_arg = _val if _is_quoted(_val, '"') else f'"{_val}"'
            tmp = container.call(
                "echo", _echo_arg, log=False, capture_output=True, **call_kwargs
            )
            if tmp.returncode != 0 or tmp.stdout[-1] != "\n":
                raise RuntimeError(
                    f"Error while executing `echo {_echo_arg}`. The returncode "
                    f"is {tmp.returncode} and output is {tmp.stdout!r}"
                )
            overrides[_var] = tmp.stdout[:-1]  # <- strip off trailing newline
        args = args[1:]

    if len(overrides):
        if call_kwargs.get("env") is None:
            # remember "env" kwarg might have been explicitly set to None
            call_kwargs["env"] = overrides
        else:
            call_kwargs["env"] = call_kwargs["env"].copy()
            call_kwargs["env"].update(overrides)

    print_info(f"$ {_msg}", flush=True)
    rslt = container.call(*args, log=False, capture_output=False, **call_kwargs)
    if expect_success and rslt.returncode != 0:
        return TestResult.Fail(_msg)
    return TestResult.Success()


def run_test(
    test_item: TestItem,
    container: DockerContainer,
    print_info: Callable = print,
) -> TestResult:
    """
    Actually executes the test

    Notes
    -----

    The Plan:
    - the basic idea here is to use:
      * ``print_info`` to describe information (e.g. useful metadata or describing a
        phase of a text) or describe an action (e.g. executing a shell command or
        copying a file)
      * ``print`` and direct usage of ``sys.stdout`` are used to communicate the result
        of an action (namely, output from executing a subprocess)
    - By default we want to be able to log all this information to a temporary
      log file while executing a test, and then only show it to the user if a
      test failure occurs
    - We separate these modes of outputs to leave an easy path to adding a
      verbosity setting where we show print_info messages while
      running and only replay the detailed log during errors

    In practice, this approach isn't optimal. We can *DEFINITELY* do better
    than this, but this is a simple enough way to start...
    """
    _cmd = functools.partial(
        _exec_container_cmd, container=container, print_info=print_info
    )

    # define the path to the root project (within the container)
    sample_proj_root = "./sample-project"

    # do the setup
    print_info(
        f"Setup: Create sample project in the container: {sample_proj_root}", flush=True
    )
    tmp = _cmd("mkdir", sample_proj_root, expect_success=True)
    if not tmp.is_success():
        return tmp
    for partial_src_path in test_item.src_files:
        src_path = os.path.join(GRACKLE_ROOT, partial_src_path)
        dst = os.path.join(sample_proj_root, os.path.basename(partial_src_path))
        pretty_src_path = os.path.join("GRACKLE", partial_src_path)
        print_info(f"-> copy {pretty_src_path}", flush=True)
        try:
            container.copy_into(src_path, dst)
        except ScriptError as e:
            return TestResult.Fail(f"copy {pretty_src_path} to CONTAINER @ {dst}", e)

    # Move onto running commands
    # ==========================

    # Execute Build Commands
    print_info("Execute Build Commands:", flush=True)
    for i, cmd in enumerate(test_item.cmds):
        is_last_build_cmd = len(test_item.cmds) - 1 == i
        expect_success = is_last_build_cmd and test_item.expect_success

        args = shlex.split(cmd)

        # actually call the command
        tmp = _cmd(*args, cwd=sample_proj_root, env=None, expect_success=expect_success)
        if not tmp.is_success():
            return tmp

    if not test_item.expect_success:
        return TestResult.Success()
    if test_item.expect_success:
        print_info("Test the resulting binary:", flush=True)
        output_binary = (
            f"{container.home_dir}/{sample_proj_root}/{test_item.output_binary}"
        )

        print_info("-> checking if the binary exists", flush=True)
        tmp = _cmd("test", "-f", output_binary, expect_success=True)
        if not tmp.is_success():
            return tmp

        if test_item.override_LD_LIBRARY_PATH:
            cur_env = container.get_environment()
            if "LD_LIBRARY_PATH" in cur_env:
                new_val = f"{_LOCAL_LIB}:{cur_env['LD_LIBRARY_PATH']}"
            else:
                new_val = _LOCAL_LIB
            print_info(f"Infer env override: LD_LIBRARY_PATH={new_val}")
            env_override = {"LD_LIBRARY_PATH": new_val}
        else:
            env_override = None

        # this is dumb...
        # -> we need to make sure that the hard-coded relative path to the data-files
        #    in compiled source file is correct
        # -> this is true if we execute the binary from src/example
        grackle_dir = test_item.image.get_grackle_dir()
        exec_dir = f"{grackle_dir}/src/example"

        print_info("-> attempting to execute the binary", flush=True)
        # if we wanted to check correctness, we could make use of the script at
        # f"{grackle_dir}/tests/scripts/code_example_checker.py"
        return _cmd(
            output_binary,
            cwd=exec_dir,
            env=env_override,
            include_outer_env=True,
            expect_success=True,
        )


class TestNameFilter:
    """
    Implements the googletest test filtering strategy

    For more details, see:
        https://google.github.io/googletest/advanced.html#running-a-subset-of-the-tests
    """

    def __init__(self, s: Optional[str] = None):
        pos_part, hyphen, neg_part = ("*" if s is None else s).partition("-")
        if s == "":
            raise ValueError("a filter string can't be an empty string")
        elif (hyphen == "-") and (neg_part == ""):
            raise ValueError("a filter string shouldn't end with a '-'")
        elif "-" in neg_part:
            raise ValueError("a filter string should have at most 1 '-'")
        self._positive_patterns = ["*"] if pos_part == "" else pos_part.split(":")
        self._negative_patterns = [] if neg_part == "" else neg_part.split(":")

    def is_match(self, name: str) -> bool:
        match_pos = any(fnmatch.fnmatch(name, pat=p) for p in self._positive_patterns)
        match_neg = any(fnmatch.fnmatch(name, pat=p) for p in self._negative_patterns)
        return match_pos and not match_neg


def collect_tests(
    scan_dir: os.PathLike, test_filter: Optional[TestNameFilter] = None
) -> Dict[str, TestItem]:
    """
    Construct a dict holding all test items
    """
    # we can be much more clever than this!
    tests = {}
    with os.scandir(scan_dir) as outer_it:
        for testdir_entry in outer_it:
            if not testdir_entry.is_dir():
                continue
            elif "-" in testdir_entry.name or ":" in testdir_entry.name:
                # this is based on the filtering logic
                raise RuntimeError("- and : aren't allowed in test case names")
            conf_path = os.path.join(testdir_entry.path, "conf.toml")
            if not os.path.isfile(conf_path):
                continue
            for name, item in test_cases_from_conf(conf_path):
                if test_filter.is_match(name):
                    tests[name] = item
    return tests
