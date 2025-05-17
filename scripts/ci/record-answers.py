#!/usr/bin/env python3
"""
Create a temporary isolated build of the specified Pygrackle revision (without
modifying the source directory) & use it to store the results of answer-tests

To achieve this, this tool creates a temporary directory ``$<basetemp>`` (with
a fresh venv inside of it) and follows the following procedure:

  1. Duplicate the Grackle repository inside of ``basetemp``
  2. Checkout the specified revision of the duplicated repository
  3. Invoke ``pip install`` with the venv and the duplicated repository
  4. Invoke ``pytest`` using the venv and duplicated repository

By default, $<basetemp> directory is cleaned up at the end
"""  # this docstring is reused for argparse's description argument

# for portability: only use standard-library modules present in older python versions
import argparse
import contextlib
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
from typing import Any, Mapping, Sequence, NamedTuple, Optional, Union

# Handle some global stuff
# ========================
if sys.version_info < (3, 6, 1):  # 3.6.0 doesn't support all NamedTuple features
    raise RuntimeError("python 3.6.1 or newer is required")

logger = logging.getLogger("recorder")
logger.setLevel(logging.DEBUG)


def _configure_logger(color=False):
    global logger
    fmt = "%(name)s > %(message)s"
    if color:
        fmt = "\x1b[36;20m" + fmt + "\x1b[0m"

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(fmt))
    logger.addHandler(console_handler)


class RecorderError(RuntimeError):
    pass


class Config(NamedTuple):
    answer_dir: str  # path to the directory where answers are recorded
    git_rev: Optional[str]  # revision to use (None corresponds to HEAD)
    color: bool  # whether to color the log messages
    interpreter: str  # python interpreter used for venv
    extra_pip_install_args: Sequence[str]  # extra args passed to pip-install
    extra_pytest_args: Sequence[str]  # extra args passed to pytest
    basetemp: Optional[str]  # explicit path to base temporary directory
    basetemp_retain: Optional[str]  # retention-policy for base temporary dir
    ref_grackle_dir: str  # path to primary grackle directory

    @classmethod
    def from_args(cls, parser, parse_args_kw=None):
        parse_args_kw = {} if parse_args_kw is None else parse_args_kw
        kw = vars(parser.parse_args(**parse_args_kw))
        for key in filter(lambda key: key not in Config._fields, kw):
            raise RuntimeError(f"{key!r} doesn't obviously map to Config")

        ref_grackle_dir = os.path.abspath(f"{os.path.dirname(__file__)}/../..")
        assert os.path.isdir(os.path.join(ref_grackle_dir, ".git"))
        out = cls(ref_grackle_dir=ref_grackle_dir, **kw)

        assert out.interpreter is not None  # sanity-check!
        return out


def _expand(arg: Optional[Union[str, Sequence[str]]], **kwds: Any) -> Any:
    # expand placeholder in arg (denoted by $<placeholder>)
    # - arg is nominally a single string or a sequence of strings.
    # - Keyword arguments specify placeholders & their replacements

    pattern = r"\$<(\w+)>"

    def repl(m: re.Match):
        try:
            return kwds[m.group(1)]
        except KeyError:
            msg = f"unknown placeholder, `$<{m.group(1)}>`, in {arg!r}"
            raise KeyError(msg) from None

    if arg is None:
        return None
    elif isinstance(arg, str):
        return re.sub(pattern, repl, arg)
    return [re.sub(pattern, repl, elem) for elem in arg]


def _run(
    *args: str,
    substitutions: Optional[Mapping[str, str]] = None,
    log: bool = True,
    silent: bool = False,
    check_returncode: bool = True,
    cwd: Optional[str] = None,
    timeout: Optional[float] = None,
    env: Optional[Mapping[str, str]] = None,
) -> int:
    """Invoke a command

    Parameters
    ----------
    *args
        The command and its arguments
    substitutions : dict, option
        Map specifying substitutions
    log : bool
        When true, we log the command
    silent : bool
        When True, stdout and stderr are suppressed
    check_returncode : bool
        When True, we report an error for non-0 return codes
    cwd : str, optional
        Optionally specifies a directory to invoke the command from
    timeout : float, optional
        If the timeout expires, the subprocess will be killed and after it
        is done terminating, an exception is raised
    env : dict, optional
        When specified, it's used to all of the subprocess's env variables
        (rather than inheriting the current process's env)

    Returns
    -------
    int
        The command's returncode
    """
    # some argument checking:
    substitutions = {} if substitutions is None else substitutions
    if len(args) == 0:
        raise ValueError("args was not specified")
    elif not isinstance(args[0], str):
        raise TypeError(f"args[0], {args[0]!r}, isn't a str")

    coerced_args = _expand(args, **substitutions)

    if log:
        _msg = " ".join(args)
        if cwd is not None:
            _msg += f"; (exec_dir: {cwd})"
        logger.info(_msg)
    rslt = subprocess.run(
        coerced_args,
        cwd=_expand(cwd, **substitutions),
        stdout=subprocess.PIPE if silent else None,
        stderr=subprocess.STDOUT if silent else None,
    )
    if check_returncode and (rslt.returncode != 0):
        if silent and rslt.stdout:
            print(rslt.stdout.decode("utf8"), file=sys.stderr, flush=True)
        cwd = "./" if cwd is None else cwd
        env = "<inherit>" if env is None else env
        raise RecorderError(
            "subprocess exited with nonzero code\n"
            f"  command: {' '.join(args)}\n  exec_dir: {cwd!r}\n"
            f"  env: {env!r}\n  code: {rslt.returncode}\n"
        )
    else:
        return rslt.returncode


def _validate_git_rev(rev, git_dir):
    # this does a quick check
    # - we could probably do a more rigorous check if we knew whether rev was
    #   supposed to be a tag/branch-name/commit-hash and if we consider some
    #   details about git-plumbing
    # - the current check is probably less restrictive than git-checkout (its
    #   important that we aren't more restrictive)

    cmd = ("git", "rev-parse", "--silent", "--verify", "--end-of-options", rev)
    returncode = _run(*cmd, silent=True, cwd=git_dir, log=False, check_returncode=False)
    if returncode == 0:
        return True
    raise RecorderError(
        f"`{rev}` is not a known git revision. If it is a tag, you may need "
        "to fetch tags from the remote repository"
    )


def _make_empty_dir(path, exist_ok=True):
    # returns True if we actually made a directory
    if not os.path.exists(path):
        os.mkdir(path)
        return True
    elif not exist_ok:
        raise RecorderError(f"{path} already exists")
    elif len(os.listdir(path)):  # may raise NotADirectoryError
        raise RecorderError(f"`{path}` isn't an empty dir")
    return False  # return False since it already existed


class _Session(NamedTuple):  # loosely inspired by nox
    """Provides access to the basetemp dir and the virtual env"""

    bin: str
    basetemp: str

    @property
    def substitutions(self):
        return {"py": self.bin, "basetemp": self.basetemp}

    def expand(self, arg):
        return _expand(arg, **self.substitutions)

    def run(self, *args, **kwargs):
        return _run(*args, substitutions=self.substitutions, **kwargs)


@contextlib.contextmanager
def open_session(conf: Config):
    # First, ensure basetemp directory exists (& ensure we have a valid policy)
    retention = conf.basetemp_retain
    assert (retention is None) or retention in ("always", "failed", "never")
    if conf.basetemp is None:
        basetemp, dflt_policy = tempfile.mkdtemp(), "never"
    else:
        path_is_newly_created = _make_empty_dir(conf.basetemp, exist_ok=True)
        if path_is_newly_created:
            basetemp, dflt_policy = conf.basetemp, "never"
        elif (retention is None) or (retention == "always"):
            basetemp, dflt_policy = conf.basetemp, "always"
        else:
            raise RecorderError(f"retention must be 'always' if `{basetemp}` exists")
    policy = retention if retention is not None else dflt_policy
    logger.info(f"using $<basetemp>: {basetemp}")

    # Next, setup the virtual environment
    _venv_dir, venv_bin = "$<basetemp>/venv", "$<basetemp>/venv/bin/python"
    logger.info(f"setup venv ($<py>: {venv_bin})")
    _cmd = [conf.interpreter, "-m", "venv", _venv_dir]
    _run(*_cmd, substitutions={"basetemp": basetemp}, log=False, silent=True)

    # Now, yield the _Session instance (and deal with cleanup)
    is_successful_exit = False
    try:
        yield _Session(bin=_expand(venv_bin, basetemp=basetemp), basetemp=basetemp)
        is_successful_exit = True
    finally:
        if (is_successful_exit and (policy == "failed")) or (policy == "never"):
            shutil.rmtree(basetemp)


def setup_and_record(session: _Session, conf: Config):
    # upgrade dependencies (may not be necessary
    _cmd = ["$<py>", "-m", "pip", "install", "--upgrade", "pip", "setuptools"]
    session.run(*_cmd, log=False, silent=True)

    if conf.git_rev is None:
        logger.info("HEAD of current grackle dir will be used")
        grackle_dir = conf.ref_grackle_dir
        if not os.path.isabs(conf.ref_grackle_dir):
            grackle_dir = f"./{grackle_dir}"
    elif True:
        # rather than cloning, I suspect it makes more sense to invoke
        # shutil.copytree and git-clean
        grackle_dir = "$<basetemp>/grackle"
        logger.info(f"cloning directory into {grackle_dir}")
        _cmd = [
            "git", "clone", "--recursive", "--branch", "main", "--tags", conf.ref_grackle_dir, grackle_dir
        ]  # fmt: skip
        session.run(*_cmd, log=True, silent=False)
        logger.info(f"checkout {conf.git_rev} from cloned directory")
        session.run(
            "git", "checkout", conf.git_rev, cwd=grackle_dir, log=True, silent=True
        )

    logger.info("install pygrackle")
    _cmd = (
        ["$<py>", "-m", "pip", "install"]
        + conf.extra_pip_install_args
        + ["-e", f"{grackle_dir}[dev]"]
    )
    session.run(*_cmd, silent=True)

    logger.info(f"record results of the gold-standard tests to {conf.answer_dir}")
    _adir = os.path.abspath(conf.answer_dir)
    test_flags = [f"--answer-dir={_adir}", "--answer-store"] + conf.extra_pytest_args
    session.run("$<py>", "-m", "pytest", *test_flags, cwd=grackle_dir)


def main(conf):
    _configure_logger(color=conf.color)

    try:
        # before anything else check on git revision & create answer-dir (to fail fast)
        if conf.git_rev is not None:
            _validate_git_rev(rev=conf.git_rev, git_dir=conf.ref_grackle_dir)
        _make_empty_dir(conf.answer_dir)

        with open_session(conf) as session:
            setup_and_record(session, conf)

    except RecorderError as err:
        logger.error(f"ERROR: {err.args[0]}")
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    except BaseException:
        logger.error("Unexpected error:")
        traceback.print_exc(file=sys.stderr)
        return 70  # https://www.man7.org/linux/man-pages/man3/sysexits.h.3head.html
    else:
        logger.info("success")
        return 0


# Define the argument parser
# ==========================
def _add_argforward_opt(parser, short, long, cmd, forbidden=None):
    # defines options for argument forwarding
    forbidden = [] if forbidden is None else forbidden

    assert long.endswith("-arg")
    _dest = f"extra_{long[2:-4].replace('-', '_')}_args"

    _help = f"Delegates args/options to `{cmd}`."
    if forbidden:
        _help = f"{_help} Forbidden args: {', '.join(forbidden)}"

    def forwarded_arg(arg):
        if any(re.match(f"^{f}($| |=)", arg) for f in forbidden):
            raise ValueError("forwarded argument is forbidden")
        elif not arg.startswith("-"):
            # this is being done so we have the option to forward ENV variables
            raise ValueError("forwarded arguments should start with -")
        return arg

    _type = forwarded_arg

    parser.add_argument(
        short, long, action="append", type=_type, dest=_dest, help=_help, default=[]
    )


parser = argparse.ArgumentParser(
    description=__doc__.strip(),  # remove leading and trailing string newlines
    formatter_class=argparse.RawDescriptionHelpFormatter,
    allow_abbrev=False,
)

parser.add_argument("answer_dir", help="path to the output answer-directory")

vergrp = parser.add_mutually_exclusive_group(required=True)
vergrp.add_argument("--git-rev", help="tag/branch/commit used to build pygrackle")
vergrp.add_argument(
    "--use-HEAD",
    dest="git_rev",
    action="store_const",
    const=None,
    help="build pygrackle from current HEAD",
)

parser.add_argument("--basetemp", help="override base temporary dir (for debug)")
parser.add_argument(
    "--basetemp-retain",
    choices=("always", "failed", "never"),
    help="Policy for retaining base temporary dir (for debug)",
)
parser.add_argument(
    "--interpreter", default=sys.executable, help="override python interpreter for venv"
)
parser.add_argument("--color", action="store_true", help="use color")

fwdgrp = parser.add_argument_group(
    "option-fowarding",
    description="""\
Each flag-set delegates options to a subprocess. We highlight scenarios
demonstrating how to use the (fictional) flags `-M,--mycmd-arg` to
delegate/forward options to invocations of the (fictional) `mycmd` command:
    (no flags)              ->   mycmd <ordinary-args...>
    -M-v                    ->   mycmd -v <oridinary-args...>
    --mycmd-arg=--foo=bar   ->   mycmd --foo=bar <ordinary-args...>
    -M-v -M--foo=bar        ->   mycmd -v --foo=bar <ordinary-args...>
Do **NOT** put a space before a delegated option starting with `-`.""",
)
_add_argforward_opt(fwdgrp, "-I", "--pip-install-arg", "pip install")
_add_argforward_opt(
    fwdgrp, "-T", "--pytest-arg", "pytest", forbidden=["--answer-dir", "--answer-store"]
)

if __name__ == "__main__":
    sys.exit(main(Config.from_args(parser)))
