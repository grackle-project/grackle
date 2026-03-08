"""
Define some useful tools used by other scripts
"""

# for portability, we restrict ourselves to built-in modules provided with python 3.6
import contextlib
import logging
import os
import re
import subprocess
import sys
from typing import Container, Dict, IO, Iterator, Mapping, NamedTuple, Optional, Union

if sys.version_info < (3, 6, 1):  # 3.6.0 doesn't support all NamedTuple features
    raise RuntimeError("python 3.6.1 or newer is required")

logger = logging.getLogger("ci-setup")
logger.setLevel(logging.DEBUG)


def configure_logger(color: bool = False):
    """
    Configure the logger for pretty outputs.

    A script should call this once, and only once.
    """
    global logger

    color_start = ""
    color_stop = ""
    if color:
        color_start = "\x1b[36;20m"
        color_stop = "\x1b[0m"

    fmt = f"{color_start}%(name)s{color_stop} > %(message)s"

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(fmt))
    logger.addHandler(console_handler)


class ScriptError(RuntimeError):
    """
    Generic error type used for expected failure modes

    Ideally, the error message should give a high-level description that
    provides the user with enough context to understand the issue without
    a traceback.
    """

    pass


if sys.version_info < (3, 10):

    def _os_release_lines() -> Iterator[str]:
        for path in ("/etc/os-release", "/usr/lib/os-release"):
            with contextlib.suppress(FileNotFoundError):
                with open(path, "r") as f:
                    yield from f
                return  # <- exit immediately
        raise OSError("unable to find '/etc/os-release' or '/usr/lib/os-release'")

    def freedesktop_os_release() -> Dict[str, str]:
        # crude backport of platform.freedesktop_os_release (which is based on
        # https://www.freedesktop.org/software/systemd/man/latest/os-release.html)

        escape_characters = [re.escape("\\"), "'", '"', "`", re.escape("$")]
        escape_pattern = re.escape("\\") + "[" + "".join(escape_characters) + "]"

        def sanitize_val_str(val):
            if (val[0] == val[-1]) and val[0] in ("'", '"'):
                val = val[1:-1]
            return re.sub(escape_pattern, lambda m: m.group(0)[1], val)

        out = {}
        for line in _os_release_lines():
            match = re.match(r"^([a-zA-Z0-9_]+)=(.*)$", line.rstrip())
            if match is not None:
                out[match.group(1)] = sanitize_val_str(match.group(2))
        return out

else:
    from platform import freedesktop_os_release


class LinuxDistroInfo(NamedTuple):
    id: str
    version_id: Optional[str]
    codename: Optional[str]

    @classmethod
    def infer_from_system(cls):
        # this will fail for non-linux systems or old/weird linux distributions
        data = freedesktop_os_release()
        id = data["ID"]  # <- standard guarantees its always defined (and lowercase)
        codename = data.get("CODENAME")
        if codename is None and id == "ubuntu":
            codename = data.get("UBUNTU_CODENAME")
        return cls(id=id, version_id=data.get("VERSION_ID"), codename=codename)


def _fmt_env_args(
    include_outer_env: bool = True,
    env: Optional[Mapping[str, str]] = None,
) -> str:
    """
    Format a string representation conveying env variables as concisely as possible
    """
    # this assumes that the env-overwrites are short
    kv_pairs = [] if env is None else (f"{k}={v}" for k, v in env.items())
    if include_outer_env and env is None:
        return "<inherit>"
    elif include_outer_env:
        return f"<inherit>.update({'; '.join(kv_pairs)})"
    elif env is None:
        return "<no-env-vars>"
    else:
        return f"{{{'; '.join(kv_pairs)}}}"


def _get_subprocess_run_env_kwarg(
    include_outer_env: bool = True,
    env: Optional[Mapping[str, str]] = None,
) -> Optional[Dict[str, str]]:
    """Construct the env kwarg for subprocess.run"""
    if include_outer_env and env is None:
        return None  # subprocess simply inherits the environment variables
    elif include_outer_env:
        out = os.environ.copy()
        out.update(env)
        return out
    elif env is None:
        return {}  # subprocess is run with no environment variables
    else:
        return env


class CmdRslt(NamedTuple):
    returncode: int  # the exit code
    stdout: Optional[str]  # the stdout stream (if captured)


def exec_cmd(
    *args: str,
    log: bool = True,
    silent: bool = False,
    cwd: Optional[str] = None,
    timeout: Optional[float] = None,
    include_outer_env: bool = True,
    env: Optional[Mapping[str, str]] = None,
    stdout: Union[IO[str], int, None] = None,
    stderr: Union[IO[str], int, None] = None,
    success_codes: Optional[Container[int]] = (0,),
    dry_run: bool = False,
) -> CmdRslt:
    """Invoke a command

    The interface is loosely inspired by the nox API

    Parameters
    ----------
    *args
        The command and its arguments
    log : bool, optional
        When True (the default), log the command.
    silent : bool
        Default is ``False``. When ``True``, silences command output and
        returns the output from this function. This is accomplished by
        combining stdout & stderr into a single stream.
    cwd : str, optional
        Optionally specifies a directory to invoke the command from
    timeout : float, optional
        If the timeout expires, the subprocess will be killed and after it
        is done terminating, an exception is raised
    include_outer_env: bool
        When True (the default), the subprocess inherits the environment of
        the current process
    env : dict, optional
        When specified, it's used to specify the subprocess's env variables.
        When include_outer_env is True, we overwrite variables.
    stdout
        Optionally specifies an open file object or a file descriptor where
        the contents of stdout are written. Incompatible with silent=True
    stderr
        Optionally specifies an open file object or a file descriptor where
        the contents of stderr are written. Incompatible with silent=True

    Returns
    -------
    CmdRslt
        Holds the return code and stdout (if it was captured)
    """
    # some argument checking:
    if len(args) == 0:
        raise ValueError("args was not specified")
    elif not isinstance(args[0], str):
        raise TypeError(f"args[0], {args[0]!r}, isn't a str")

    if log:
        _msg = " ".join(args)
        _meta_list = []
        if cwd is not None:
            _meta_list.append(f"exec_dir: {cwd}")
        _env_str = _fmt_env_args(include_outer_env=include_outer_env, env=env)
        _meta_list.append(f"ENV: {_env_str}")
        logger.info(f"$ {_msg}; ({'; '.join(_meta_list)})")

    if dry_run:
        return

    # adjust stdout & stderr if necessary
    if silent:
        if stdout is not None:
            raise ValueError("Can't specify stdout kwarg with silent==True")
        elif stderr is not None:
            raise ValueError("Can't specify stderr kwarg with silent==True")
        # combine stdout and stder into a single stream
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT
    elif stderr == subprocess.PIPE:
        raise ValueError("currently no support for stderr=subprocess.PIPE")

    tmp_rslt = subprocess.run(
        args,
        cwd=cwd,
        stdout=stdout,
        stderr=stderr,
        env=_get_subprocess_run_env_kwarg(include_outer_env=include_outer_env, env=env),
        timeout=timeout,
    )
    sys.stdout.flush()

    # repackage the result
    _stdout = tmp_rslt.stdout.decode("utf8") if tmp_rslt.stdout is not None else None
    rslt = CmdRslt(returncode=tmp_rslt.returncode, stdout=_stdout)

    if (success_codes is not None) and (rslt.returncode not in success_codes):
        if silent and rslt.stdout:
            print(rslt.stdout, file=sys.stderr, flush=True)
        cwd = "./" if cwd is None else cwd
        raise ScriptError(
            "subprocess exited with nonzero code\n"
            f"  command: {' '.join(args)}\n  exec_dir: {cwd!r}\n"
            f"  env: {_fmt_env_args(include_outer_env=include_outer_env, env=env)}\n"
            f"  code: {rslt.returncode}\n"
        )
    return rslt
