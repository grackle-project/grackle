"""
Define the core logic for running install tests
"""

# to maximize portability, only use built-in modules present in python 3.7 or newer
from enum import Enum
import fnmatch
from functools import partial
import os
import re
import shlex
import textwrap
from types import MappingProxyType
from typing import Callable, Dict, Iterator, Mapping, NamedTuple, Optional, Tuple, Union

# import local logic
import conf_reader
from docker import DockerContainer, GRACKLE_ROOT, import_standalone_module


common = import_standalone_module(
    dir_path=os.path.join(GRACKLE_ROOT, "scripts", "ci"), module_name="common"
)
ScriptError = common.ScriptError
logger = common.logger
make_cmd_summary_str = common.make_cmd_summary_str

DOCKERFILE_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "installtest.Dockerfile"
)

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


def _build_req_mapping(
    *reqs: conf_reader.TOMLValReq,
) -> Mapping[str, conf_reader.TOMLValReq]:
    out = {req.name: req for req in reqs}
    assert len(out) == len(reqs)
    return MappingProxyType(out)


_SHARED_TAB_REQS = _build_req_mapping(
    conf_reader.TOMLValReq(
        name="src_file",
        classinfo=str,
        type_descr="string",
        description=textwrap.dedent(
            """\
           the file from GRACKLE_ROOT/src/example to copy into the sample
           project's source directory."""
        ),
    ),
    conf_reader.TOMLValReq(
        name="output_binary",
        classinfo=str,
        type_descr="string",
        required=False,
        description=textwrap.dedent(
            """\
           The relative path to the output binary that we expect to be created
           from the specified commands. This path is relative to the root of
           the sample project's source directory."""
        ),
    ),
)

_TEST_TAB_REQS = _build_req_mapping(
    conf_reader.PARAMETRIZED_TOML_REQ,
    conf_reader.TOMLValReq(
        name="image",
        classinfo=str,
        type_descr="string",
        choices=set(member.layer_name() for member in DockerImage),
        description=textwrap.dedent(
            """\
            Specifies the names of the images that is used to initialize the
            container in which the test is run. Essentially, you can think of
            the image as a cached set of setup steps."""
        ),
    ),
    conf_reader.TOMLValReq(
        name="cmds",
        classinfo=partial(
            conf_reader._is_single_type_toml_arr, kind=str, allow_empty=True
        ),
        type_descr="array of strings",
        description=textwrap.dedent(
            """\
            Each entry of the array specifies a separate command that is
            executed in a sequence.

            Each command is executed in a separate subshell. Practically, this
            means that modifying an environment variable in one command won't
            affect the value during subsequent commands

            .. note::

               In the future, it may be worth allowing array to contain a mixture
               of strings and inline tables. The could make it easier for us to
               handle environment variables. We illustrated the difference below:

               .. code-block:: toml

                  # current approach
                  cmds = [
                    "Grackle_DIR='/home/gr-user/local' cmake -Bbuild",
                    "cmake --build build"
                  ]
                  # alternative approach approach
                  cmds = [
                    {cmd = "cmake -Bbuild", env.Grackle_DIR="/home/gr-user/local"},
                    "cmake --build build"
                  ]

               We probably don't want to do this unless we support TOML 1.1 (since TOML
               1.0 requires inline-tables to be fully defined on a single line)."""
        ),
    ),
    conf_reader.TOMLValReq(
        name="expect_cmd_fail",
        classinfo=bool,
        type_descr="boolean",
        required=False,
        description=textwrap.dedent(
            """\
            Indicates whether we expect the final command to fail. Omitting
            this parameter is equivalent to setting it to ``false``. When
            ``true``, the output_binary parameter can be omitted"""
        ),
    ),
    conf_reader.TOMLValReq(
        name="override_LD_LIBRARY_PATH",
        classinfo=bool,
        type_descr="boolean",
        required=False,
        description=textwrap.dedent(
            """\
            Specify whether the ``LD_LIBRARY_PATH`` environment variable needs
            to be overriden in order to execute the output binary. Omitting
            this parameter is equivalent to setting it to ``false``

            When performing building a test program manually (e.g. with a
            Makefiles), this is typically required. This generally isn't needed
            when a test program with CMake unless it's both linked against a
            shared library AND the program is installed (if executing the test
            program from the build-directory, this isn't needed)."""
        ),
    ),
)


def test_cases_from_conf(file_path: os.PathLike) -> Iterator[Tuple[str, TestItem]]:
    """Constructs TestItem from the specified configuration toml file"""

    # before we actually parse the config file, let's build a list of all other
    # files in the directory that holds the config file
    # -> we assume that these are parts of the sample source directory
    # -> we store the paths to these files so that they're relative to GRACKLE_ROOT
    abs_conf_dir_path = os.path.dirname(os.path.abspath(file_path))
    _rel_conf_dir_path = os.path.relpath(abs_conf_dir_path, start=GRACKLE_ROOT)
    common_src_files = []
    with os.scandir(abs_conf_dir_path) as it:
        for entry in it:
            if entry.is_dir(follow_symlinks=True):
                raise ValueError(f"dir holding {file_path} can't hold subdirectories")
            elif os.path.samefile(entry.path, file_path):
                continue  # <- we won't need to copy the conf file itself
            else:
                common_src_files.append(os.path.join(_rel_conf_dir_path, entry.name))
    # later on, we will use conf_dirname as the prefix for each test name
    conf_dirname = str(os.path.basename(abs_conf_dir_path))

    reader = conf_reader.ConfReader(file_path)

    if len(reader.toml_doc) > 2:
        extra = list(filter(lambda k: k != "shared" and k != "test", reader.toml_doc))
        raise conf_reader.ConfReadError(
            f"unexpected table(s): {', '.join(extra)}", reader
        )

    # read the shared table
    shared_table, _ = reader.process_table(path="shared", reqs=_SHARED_TAB_REQS)
    src_file_path = os.path.join("src", "example", shared_table["src_file"])
    src_files = tuple([src_file_path] + common_src_files)

    test_table_paths = [("test", name) for name in reader["test"].keys()]
    for table_path in test_table_paths:
        # read the test table
        test_table, parametrization = reader.process_table(table_path, _TEST_TAB_REQS)

        # construct each variation
        n_variations = 1 if parametrization is None else len(parametrization)
        for i in range(n_variations):
            if parametrization is None:
                test_name = f"{conf_dirname}.{table_path[-1]}"
                _table = test_table
            else:
                test_name = f"{conf_dirname}.{table_path[-1]}/{i}"
                _table = {}
                for k, v in test_table.items():
                    if isinstance(v, conf_reader.ParametrizedProxy):
                        _table[k] = v.substitute(reader, parametrization[i])
                    else:
                        _table[k] = v

            image = DockerImage(_table["image"])

            if _table.get("expect_cmd_fail", False):
                output_binary = None
            elif "output_binary" not in shared_table:
                msg = "when True (or not specified), shared.output_binary is required"
                full_path = table_path + ("expect_cmd_fail",)
                raise conf_reader.ConfReadError(msg, reader, path=full_path)
            else:
                output_binary = shared_table["output_binary"]

            test_item = TestItem(
                image=image,
                src_files=src_files,
                cmds=_table["cmds"],
                output_binary=output_binary,
                override_LD_LIBRARY_PATH=_table["override_LD_LIBRARY_PATH"],
            )
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
