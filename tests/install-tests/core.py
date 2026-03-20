"""
Define the core logic for running install tests
"""

# to maximize portability, only use built-in modules present in python 3.6.1 and newer
from collections import ChainMap
from enum import Enum
import fnmatch
from functools import partial
import os
import re
import shlex
import sys
from typing import (
    Any,
    Callable,
    Container,
    Dict,
    Iterator,
    List,
    Mapping,
    NamedTuple,
    Optional,
    Tuple,
    Type,
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
        return f"gr_installtest_{self._value_}"

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

_CheckTypeArg = Union[Type, Tuple[Type, ...], Callable[[Any], bool]]


def _check_type(v: Any, classinfo: _CheckTypeArg) -> bool:
    if isinstance(classinfo, (type, tuple)):
        return isinstance(v, classinfo)
    return classinfo(v)


class TomlValReq(NamedTuple):
    """
    Specifies requirements for a toml value

    Checks are applied after the values are coerced to python values.
    """

    # this is either an argument for isinstance or a callable specifying whether
    # the value has the appropriate type
    classinfo: _CheckTypeArg

    # describes the required toml type
    type_descr: str

    # when specified, the value is required to be one of these values
    choices: Optional[Container] = None

    # when True, the associated parameter is allowed as an array
    # -> this is intended as a quick-and-dirty mechanism for parameterizing
    #    several related tests
    # -> as soon as we have 2 or more parameters where we want to allow
    #    parameterization (if it ever happens), I **REALLY** think we want to
    #    transition to a system similar to github actions, where we introduce
    #    a matrix subtable
    #    https://docs.github.com/en/actions/how-tos/write-workflows/choose-what-workflows-do/run-job-variations
    allow_as_arr: bool = False


def _is_single_type_toml_arr(v: Any, kind: _CheckTypeArg, allow_empty: bool) -> bool:
    """Check if ``v`` corresponds to a toml array of a single type"""
    if not isinstance(v, list):
        return False
    elif len(v) == 0:
        return allow_empty
    return all(_check_type(elem, kind) for elem in v)


_is_str_arr = partial(_is_single_type_toml_arr, kind=str, allow_empty=True)

# lists the requirements for all allowed keys in one of the TOML test tables
_REQS: Dict[str, TomlValReq] = {
    # the name of the image that we want to use
    "image": TomlValReq(
        classinfo=str,
        type_descr="string",
        choices=set(member.layer_name() for member in DockerImage),
        allow_as_arr=True,
    ),
    # the list of commands that we want to use
    # -> if we tweak the vendored toml parser to support toml 1.1, I think it
    #    would be better if expected an array of (inline) tables (where we could
    #    override environment variables) or a array of mixed strings and tables
    "cmds": TomlValReq(_is_str_arr, "array of strings"),
    # the file from GRACKLE_ROOT/src/example to copy into the sample project
    # source directory
    "src_file": TomlValReq(str, "string"),
    # the relative path to the output binary (from the root of the sample project
    # source directory)
    "output_binary": TomlValReq(str, "string"),
    # indicates whether we expect the final command to fail. When this is True,
    # output_binary can be omitted
    "expect_cmd_fail": TomlValReq(bool, "boolean"),
    # specify whether we need to override_LD_LIBRARY_PATH in order to execute the
    # output binary
    "override_LD_LIBRARY_PATH": TomlValReq(bool, "boolean"),
}


def _process_toml_table(
    table: Mapping, is_common_table: bool, err_prefix=""
) -> Tuple[Dict[str, Any], Optional[Tuple[str, List[Any]]]]:
    """
    Process a (already parsed) toml table

    Parameters
    ----------
    table
        The already parsed toml table
    is_common_table
        Pass ``True`` if ``table`` corresponds to the toml file's `common`
        table. Otherwise, pass ``False``.
    err_prefix
        An optional string to prepend to error messages

    Returns
    -------
    sanitized: dict
        A sanitized dict where values are guaranteed to satisfy all
        requirements. This excludes any parametrized values
    parametrization: tuple or None
        If table included a parametrized value, this is a pair. The first
        element is the parameter name. The second is list of values.
    """

    def _mk_err(msg):
        return ValueError(err_prefix + msg)

    parametrization = None
    sanitized = {}

    for key, val in table.items():
        try:
            req = _REQS[key]
        except KeyError:
            raise _mk_err(f"unexpected parameter: {key}") from None

        if _check_type(val, req.classinfo):
            if (req.choices is not None) and (val not in req.choices):
                raise _mk_err(f"{key} must hold a value in {list(req.choices)!r}")
            sanitized[key] = val
        elif not req.allow_as_arr:
            raise _mk_err(f"{key} must have type: {req.type_descr}")
        else:
            type_is_ok = _is_single_type_toml_arr(val, req.classinfo, allow_empty=True)
            if type_is_ok and is_common_table:
                # there is a reason that this MUST be enforced
                raise _mk_err(f"common.{key} is forbidden from being an array")
            elif type_is_ok and len(val) == 0:
                raise _mk_err(f"`{key}` must not be an empty array")

            valid = (req.choices is None) or all(elem in req.choices for elem in val)
            if not valid:
                raise _mk_err(
                    f"when {key} is an array, each array element must be a value from "
                    f"{list(req.choices)!r}"
                )
            elif len(set(val)) != len(val):
                raise _mk_err(f"when {key} is an array, array element can repeat")
            elif parametrization is not None:
                # if this ever comes to pass, we probably want to rethink how things
                # work. I outline an alternative plan in the docstring comment of
                # TomlValReq's allow_as_arr attribute
                raise NotImplementedError(
                    "can't handle the case when 2 parameters are parametrized"
                )
            parametrization = (key, val)
    return sanitized, parametrization


def test_cases_from_conf(path: str) -> Iterator[Tuple[str, TestItem]]:
    """Constructs TestItem from the specified configuration toml file"""

    # before we actually parse the config file, let's build a list of all other
    # files in the directory that holds the config file
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

    with open(path, "rb") as f:
        raw = tomllib.load(f)

    err_prefix = f"Toml File Error {path}: "

    # read top-level entities and do error-checks
    common_table = raw.get("common", {})
    if not isinstance(common_table, dict):
        raise ValueError(f'{err_prefix}"common" key must provide a table, if specified')
    common_table, _dummy = _process_toml_table(
        common_table, is_common_table=True, err_prefix=err_prefix
    )
    assert _dummy is None, "sanity check!"

    test_tables = raw.get("test", {})
    if not isinstance(test_tables, dict):
        raise ValueError(f'{err_prefix}"test" is an invalid toml key')
    elif len(test_tables) == 0:
        raise ValueError(f"{err_prefix}no test tables defined")

    name_prefix = os.path.basename(abs_conf_dir_path)
    for table_name, raw_table in test_tables.items():
        if not isinstance(raw_table, dict):
            raise ValueError(f"{err_prefix}test.{table_name} is not a table")
        sanitized, parametrization = _process_toml_table(
            raw_table, is_common_table=False, err_prefix=err_prefix
        )
        base_name = f"{name_prefix}.{table_name}"
        base_conf = ChainMap(sanitized, common_table, TestItem._field_defaults)
        if parametrization is None:
            name_conf_pairs = [(base_name, base_conf)]
        else:
            p_name, p_vals = parametrization
            if p_name in common_table:
                # we're intentionally being a little strict here (we can relax later)
                raise ValueError(
                    f"since test.{table_name}.{p_name} is parametrized, it's an error "
                    f"to specify common.{p_name}"
                )
            name_conf_pairs = (
                (f"{base_name}/{v}", base_conf.new_child({p_name: v})) for v in p_vals
            )

        for test_name, test_conf in name_conf_pairs:
            for key in ["image", "src_file", "cmds"]:
                if key not in test_conf:
                    raise ValueError(
                        f'{err_prefix}"{key}" is missing from the specification of the '
                        f'"{test_name}" test'
                    )
            src_file_path = os.path.join("src", "example", test_conf["src_file"])

            kwargs = {
                k: v
                for k, v in test_conf.items()
                if k not in ["expect_cmd_fail", "src_file"]
            }
            kwargs["src_files"] = tuple([src_file_path] + common_src_files)
            kwargs["image"] = DockerImage(kwargs["image"])
            kwargs["cmds"] = tuple(kwargs["cmds"])
            if test_conf.get("expect_cmd_fail", False):
                kwargs["output_binary"] = None
            elif "output_binary" not in kwargs:
                raise ValueError(
                    f'{err_prefix}"{test_name}" test is missing output_binary '
                    "parameter (this is only allowed if expect_cmd_fail=true)"
                )

            test_item = TestItem(**kwargs)

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
    elif (not expect_success) and rslt.returncode == 0:
        return TestResult.Fail(_msg, issue="expected failure, but got returncode of 0")
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
    _cmd = partial(_exec_container_cmd, container=container, print_info=print_info)

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
        expect_success = (not is_last_build_cmd) or test_item.expect_success

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
