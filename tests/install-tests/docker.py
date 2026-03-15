"""
Define a DockerContainer context manager. When you enter the context, a docker
container is launched and the object associated with the context is used to
run commands within the container.
"""

import importlib.util
import io
import json
import os
import shlex
import subprocess
import sys
import uuid
from typing import Dict, IO, Mapping, Optional


# handle a few extra imports
# ==========================
def import_standalone_module(dir_path, module_name):
    path = os.path.join(dir_path, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


GRACKLE_ROOT = os.path.abspath(
    os.path.normpath(os.path.join(os.path.dirname(__file__), "..", ".."))
)
common = import_standalone_module(
    dir_path=os.path.join(GRACKLE_ROOT, "scripts", "ci"), module_name="common"
)

CmdRslt = common.CmdRslt
exec_cmd = common.exec_cmd
logger = common.logger
make_cmd_summary_str = common.make_cmd_summary_str


# define some helper functions that are used within DockerContainer
# =================================================================
def _call(
    *args: str,
    container_bash_process: subprocess.Popen,
    default_home_dir: str,
    log: bool = True,
    capture_output: bool = False,
    include_outer_env: bool = True,
    env: Optional[Mapping[str, str]] = None,
    cwd: Optional[str] = None,
) -> CmdRslt:
    """
    calls a program within a container

    IMPORTANT NOTES:
      1. this was copied from cibuildwheel and then modified
         https://github.com/pypa/cibuildwheel/blob/fd27532d12a7817f5e77461d90171d87f2a1092b/cibuildwheel/oci_container.py#L428
      2. this has been separated out of the ``DockerContainer`` class
         precisely because NONE of the code within DockerContainer was
         directly copied
    """
    chdir = "" if cwd is None else f"cd {cwd}"
    env_assignments = (
        " ".join(f"{shlex.quote(k)}={shlex.quote(v)}" for k, v in env.items())
        if env is not None
        else ""
    )
    if not include_outer_env:
        sep = "" if env_assignments == "" else " "
        env_assignments = f"-i{sep}{env_assignments}"

    command = " ".join(shlex.quote(str(a)) for a in args)
    end_of_message = "GrInstallTestSentinel" + str(uuid.uuid4())

    # log the command we're executing
    if log:
        _msg = make_cmd_summary_str(
            args=args,
            include_outer_env=include_outer_env,
            env=env,
            cwd=cwd,
            inherit_string="<container-inherit>",
        )
        logger.info(f"container $ {_msg}")

    # Write a command to the remote shell. First we change the
    # cwd, if that's required. Then, we use the `env` utility to run
    # `command` inside the specified environment. We use `env` because it
    # can cope with spaces and strange characters in the name or value.
    # Finally, the remote shell is told to write a footer - this will show
    # up in the output so we know when to stop reading, and will include
    # the return code of `command`.
    container_bash_process.stdin.write(
        bytes(
            f"""(
        {chdir}
        env {env_assignments} {command} 2>&1
        printf "%04d%s\n" $? {end_of_message}
    )
    """,
            encoding="utf8",
            errors="surrogateescape",
        )
    )
    container_bash_process.stdin.flush()
    if capture_output:
        output_io: IO[bytes] = io.BytesIO()
    else:
        output_io = sys.stdout.buffer

    while True:
        line = container_bash_process.stdout.readline()

        if line.endswith(bytes(end_of_message, encoding="utf8") + b"\n"):
            # fmt: off
            footer_offset = (
                len(line)
                - 1  # newline character
                - len(end_of_message)  # delimiter
                - 4  # 4 return code decimals
            )
            # fmt: on
            return_code_str = line[footer_offset : footer_offset + 4]
            return_code = int(return_code_str)
            # add the last line to output, without the footer
            output_io.write(line[0:footer_offset])
            output_io.flush()
            break
        else:
            output_io.write(line)
            output_io.flush()

    if isinstance(output_io, io.BytesIO):
        output = str(output_io.getvalue(), encoding="utf8", errors="surrogateescape")
    else:
        output = ""

    return CmdRslt(returncode=return_code, stdout=output)


def _get_environment(
    container_bash_process: subprocess.Popen, default_home_dir: str
) -> Dict[str, str]:
    """
    Gets the environment within the container

    IMPORTANT NOTES:
      1. this was copied from cibuildwheel and then modified
         https://github.com/pypa/cibuildwheel/blob/fd27532d12a7817f5e77461d90171d87f2a1092b/cibuildwheel/oci_container.py#L510
      2. this has been separated out of the ``DockerContainer`` class
         precisely because NONE of the code within DockerContainer was
         directly copied
    """
    one_liner = "import sys, json, os; json.dump(os.environ.copy(), sys.stdout)"
    _cmd = ("python3", "-c", one_liner)
    return json.loads(
        _call(
            *_cmd,
            container_bash_process=container_bash_process,
            default_home_dir=default_home_dir,
            log=False,
            capture_output=True,
        ).stdout
    )


# to be called when creating a context image
class DockerContainer:
    """
    This is a context manager for interacting with a container

    History
    -------
    There are essentially 3 steps:
    1. initially:
       - we created a container, launched a python script within the container
         the python script container ran the relevant commands, and then we
         waited for the container/script to end.
       - this worked, but it got clunky... the logic got increasingly complex
         once I updated the Dockerfile to only copy part of the grackle
         repository during image creation (and we copied the rest before
         starting the container)
         -> to be clear, this choice makes a lot of sense. Basically, we want
            to avoid copying the logic for running the install-tests since
            each change to the install-test logic would invalidate the docker
            cache (slowing iterative testing down dramatically)
       - I had not yet figured out a way to cleanly communicate the results
         of individual commands back to the test-driver...
    2. I toyed around with creating a tiny server program (in python) that
       would get launched in the container:
       - The premise is that it would receive commands over STDIN, and stream
         back the results over STDOUT. Importantly, it would break up the
         output into chunked messages that began with a tiny header that
         specified the type of message and its length
       - I never got far implementing this approach (it's the most option but
         would take the most work)
    3. The current approach (inspired by cibuildwheel).
       - this is similar in spirit to the previous approach, but MUCH simpler
       - essentially, we start the container by invoking bash and just stream
         bash commands directly to the container and directly receive the
         results
       - the only downside is that the approach for communicating returncodes
         is a hack (it's good enough for us) and I think there could be
         problems with extremely long STDOUT lines (but I think that would be
         a problem for any pipelined program -- it shouldn't be an issue)

    Implementation Note
    -------------------
    The contents of this class were written independently of the contents of
    cibuildwheel. I had already figured out the docker commands myself, and it
    was a simple matter of refactoring and slightly adapting them. This
    obviously bears a bunch of resemblence to cibuildwheel's OCIContainer class,
    but that comes from natural convergence
    """

    image_name: str
    debug_log: bool
    # the following are only initialized after we create the container
    container_proc: Optional[subprocess.Popen] = None
    container_id: Optional[str]
    home_dir: Optional[str]

    def __init__(self, image_name: str, debug_log: bool = False):
        self.image_name = image_name
        self.debug_log = debug_log
        self.container_proc = None
        self.container_id = None
        self.home_dir = None

    def __enter__(self):
        _cmd = (
            "docker",
            "container",
            "create",
            "--interactive",
            self.image_name,
            "/bin/bash",
        )
        rslt = exec_cmd(*_cmd, log=self.debug_log, silent=True)
        self.container_id = rslt.stdout.rstrip()

        _cmd = (
            "docker",
            "inspect",
            "--format={{.Config.WorkingDir}}",
            self.container_id,
        )
        self.home_dir = exec_cmd(*_cmd, log=self.debug_log, silent=True).stdout.rstrip()

        _cmd = (
            "docker",
            "container",
            "start",
            "--interactive",
            "--attach",
            self.container_id,
        )
        process = subprocess.Popen(_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        self.container_proc = process.__enter__()

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            self.container_proc.__exit__(exc_type, exc_value, traceback)
        finally:
            # let's cleanup the image
            _cmd = ("docker", "container", "rm", self.container_id)
            exec_cmd(*_cmd, log=self.debug_log, silent=True)

    def call(
        self,
        *args: str,
        log: bool = True,
        capture_output: bool = False,
        include_outer_env: bool = True,
        env: Optional[Mapping[str, str]] = None,
        cwd: Optional[str] = None,
    ) -> CmdRslt:
        return _call(
            *args,
            container_bash_process=self.container_proc,
            default_home_dir=self.home_dir,
            log=log,
            include_outer_env=include_outer_env,
            env=env,
            cwd=cwd,
        )

    def get_environment(self) -> Dict[str, str]:
        return _get_environment(
            container_bash_process=self.container_proc, default_home_dir=self.home_dir
        )

    def copy_into(self, host_path, container_path):
        if not container_path.startswith("/"):
            # container_path must be an absolute path
            container_path = os.path.join(self.home_dir, container_path)

        exec_cmd(
            "docker",
            "container",
            "cp",
            host_path,
            f"{self.container_id}:{container_path}",
            log=self.debug_log,
            silent=True,
        )
