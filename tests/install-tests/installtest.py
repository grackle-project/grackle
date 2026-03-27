#!/usr/bin/env python3
"""
Drives all of the test-cases

This will use docker to build all the relevant images and then for each test
case, it runs the run_test function
"""

import argparse
import contextlib
from datetime import datetime
import functools
import itertools
import os
import subprocess
import sys
import tempfile
import textwrap
import traceback
from typing import Any, Callable

from core import (
    DockerImage,
    TestNameFilter,
    TestResult,
    collect_tests,
    run_test,
    DOCKERFILE_PATH,
    _TOML_PARAM_REQS,
)
from docker import (
    DockerContainer,
    GRACKLE_ROOT,
    import_standalone_module,
    extract_target_info,
)

if sys.version_info < (3, 7):
    raise RuntimeError("python 3.7 or newer is required")


common = import_standalone_module(
    dir_path=os.path.join(GRACKLE_ROOT, "scripts", "ci"), module_name="common"
)
configure_logger = common.configure_logger
exec_cmd = common.exec_cmd
logger = common.logger
make_cmd_summary_str = common.make_cmd_summary_str


def _mk_docker_image(
    image: DockerImage, container: None = None, print_info: Callable = print
) -> TestResult:
    # container is just a dummy argument here

    _cmd = (
        "docker",
        "build",
        "--progress=plain",
        "--file=tests/install-tests/installtest.Dockerfile",
        f"--tag={image.tag_name()}",
        f"--target={image.layer_name()}",
        ".",
    )
    cwd = GRACKLE_ROOT
    _msg = make_cmd_summary_str(_cmd, cwd=cwd)
    print_info(f"$ {_msg}")
    tmp = exec_cmd(
        *_cmd, cwd=cwd, log=False, stderr=subprocess.STDOUT, stdout=sys.stdout
    )
    if tmp.returncode != 0:
        return TestResult.Fail(_msg, f"returncode is {tmp.returncode}")
    return TestResult.Success()


def _exec_main(args: argparse.Namespace):
    """The main function for "exec-mode", which actually drives the tests."""
    test_cases = collect_tests(
        scan_dir=args.scan_dir, test_filter=TestNameFilter(args.filter_tests)
    )

    if args.list_tests:
        for test_case in test_cases:
            print(test_case)
        return 0

    # determine the docker layers that we need
    if args.skip_image_creation:
        images = []
    else:
        images = set(case.image for case in test_cases.values())
    n_tasks = len(images) + len(test_cases)

    task_itr = itertools.chain(images, test_cases.items())

    print(f"Executing {n_tasks} tasks (including image construction)")

    # Make some preparations for redirection within tasks. Within a task, we use
    # - ``print_info`` to output information (e.g. useful metadata, describe a
    #   phase of a test) or describe an action (e.g. execute a shell command or
    #   copy a file)
    # - ``print`` and ``sys.stdout`` are used to communicate the result of an action
    #   (namely, output from executing a subprocess)
    if args.verbosity_level == 0:
        # redirected everything to a temporary buffer and only show it upon failure
        tmp_buf_cm = tempfile.TemporaryFile(mode="w+")
        print_info = print
    elif args.verbosity_level == 1:
        # like the last case we the full detailed output to a temporary buffer that
        # is revealed only upon failure. However, in this case, we also show the
        # output of print_info as the tests are run
        tmp_buf_cm = tempfile.TemporaryFile(mode="w+")

        def print_info(*args, flush: bool = True):
            # this first call gets redirected to the tmp_buf
            # (if we omit this call, then tmp_buf will be missing important context)
            print(*args, file=sys.stdout, flush=flush)
            # this second call isn't redirected
            print(*args, file=sys.__stdout__, flush=flush)

    else:
        # there's no redirection (everything is goes directly to stdout)
        tmp_buf_cm = contextlib.nullcontext()
        print_info = print

    fail_count = 0
    with tmp_buf_cm as tmp_buf:
        for task_num, task in enumerate(task_itr, start=1):
            if isinstance(task, DockerImage):
                task_name = f"make_docker_image.{task.tag_name()}"
                fn = functools.partial(_mk_docker_image, image=task)
                execution_image_name = None
            else:
                task_name, test_item = task
                fn = functools.partial(run_test, test_item=test_item)
                execution_image_name = test_item.image.tag_name()

            print(f"{task_num}/{n_tasks}: {task_name}", flush=True)

            with contextlib.ExitStack() as stack:
                if execution_image_name is None:
                    container = None
                else:
                    container = stack.enter_context(
                        DockerContainer(image_name=execution_image_name)
                    )
                if tmp_buf is not None:
                    tmp_buf.seek(0)
                    tmp_buf.truncate()
                    stack.enter_context(contextlib.redirect_stdout(tmp_buf))
                t0 = datetime.now()
                rslt = fn(container=container, print_info=print_info)
                elapsed = datetime.now() - t0
            if rslt.is_success():
                print("SUCCESS", flush=True)
                print("elapsed:", elapsed)
                continue
            print("FAILURE")
            print("-> showing the test commands and output", flush=True)
            tmp_buf.seek(0)
            for line in tmp_buf:
                print(">>", line, end="")
            print("", flush=True)

            if rslt.action_during_failure is None:
                print("Programming Error occured:\n")
                print("Exception:", rslt.issue)
                print("Traceback:")
                traceback.print_exc(file=sys.stderr)
                sys.exit(1)
            print("Error Detected During:", rslt.action_during_failure)
            if rslt.issue is not None:
                print("Issue:", rslt.issue)

            fail_count += 1

    print(80 * "=")
    print("SUMMARY:")
    print(f"{fail_count} failures")

    if fail_count == 0:
        return 0
    return 1


def _setup_exec_subcommand(subparsers: Any):
    """Configure the exec subcommand"""
    p = subparsers.add_parser("exec", help="collect & execute the full set of tests")
    p.set_defaults(fn=_exec_main)
    p.add_argument("--list-tests", action="store_true", help="list all test names")
    p.add_argument(
        "--scan-dir",
        action="store",
        default=os.path.join(GRACKLE_ROOT, "tests", "install-tests", "cases"),
        help="directory to scan for tests",
    )
    _filter_help = (
        "filters test names using the approach of googletest. For more detail, see "
        "https://google.github.io/googletest/advanced.html#running-a-subset-of-the-tests"
    )
    p.add_argument("--filter-tests", action="store", help=_filter_help)

    # TODO: figure out weirdness with creating docker images so we can remove this
    #   -> when we first build docker images, we get the expected caching of build
    #      stages that are common to all of the images
    #   -> on subsequent rebuilds, some expensive steps are not cached unless we
    #      delete all of the images first
    parser.add_argument(
        "--skip-image-creation",
        action="store_true",
        help=(
            "tells the executor to assume that all docker images were previously created "
            "and are up to date"
        ),
    )


def _entrydoc_main(args: argparse.Namespace):
    """The main function for "entrydoc-mode", which determines the"""

    packs = []  # <- we'll fill this up
    if args.kind == "param":
        for name, toml_val_req in _TOML_PARAM_REQS.items():
            extra = {"type_descr": toml_val_req.type_descr}
            if toml_val_req.choices is not None:
                extra["choices"] = ", ".join(
                    f'"{e}"' for e in sorted(toml_val_req.choices)
                )

            packs.append(
                {"name": name, "descr": toml_val_req.description, "extra": extra}
            )

    elif args.kind == "image-target":
        raw_info = extract_target_info(path=DOCKERFILE_PATH)

        # let's confirm that the DockerImage enum is exhaustive
        names_from_enum = set(e.layer_name() for e in DockerImage.__members__.values())
        diff = names_from_enum.symmetric_difference(raw_info.keys())
        if len(diff) != 0:
            raise RuntimeError(
                "The docker layer names in the DockerLayer Enum are:\n"
                f"  {', '.join(sorted(names_from_enum))}\n"
                "The docker layer names from parsing the docker file are:\n"
                f"  {', '.join(sorted(raw_info))}\n"
                "The following layer names are in set but not the other:\n"
                f"  {', '.join(sorted(diff))}"
            )
        for enumeration in DockerImage.__members__.values():
            layer_name = enumeration.layer_name()
            info = raw_info[layer_name]

            descr = textwrap.dedent("\n".join(info.description))
            extra = {
                "image-name": enumeration.tag_name(),
                # we only record dependency if it's another layer that we build
                "dependency": info.dependency if info.dependency in raw_info else "N/A",
            }
            packs.append(
                {"name": enumeration.layer_name(), "descr": descr, "extra": extra}
            )
    else:
        raise RuntimeError(f"unexpected kind argument: {args.kind!r}")

    for i, pack in enumerate(sorted(packs, key=lambda p: p["name"])):
        if i != 0:
            print("\n\n", end="")
        directive_body_indent = "   "
        pack_lines = [f".. object:: {pack['name']}", ""]
        for fieldname, value in pack["extra"].items():
            value = value.strip()
            assert "\n" not in fieldname
            assert "\n" not in value
            pack_lines.append(f"{directive_body_indent}:{fieldname}: {value.strip()}")
        pack_lines.append("")
        pack_lines.append(textwrap.indent(pack["descr"], directive_body_indent))
        print(*pack_lines, sep="\n")


def _entrydoc_subcommand(subparsers: Any):
    """Configure the entrydoc subcommand"""
    p = subparsers.add_parser(
        "entrydoc",
        help="generate rst reference documentation for each entry of the specified kind",
    )
    p.set_defaults(fn=_entrydoc_main)
    p.add_argument(
        "kind",
        choices=["param", "image-target"],
        help="the kind of entries to document",
    )


def main(args: argparse.Namespace):
    configure_logger(color=True)
    fn = args.fn
    return fn(args)


parser = argparse.ArgumentParser(
    description="drive one or more test cases",
    # prior to 3.8, setting allow_abbrev=False messed up interpretation of -vv
    allow_abbrev=sys.version_info < (3, 8),
)

subparsers = parser.add_subparsers()
_setup_exec_subcommand(subparsers)
_entrydoc_subcommand(subparsers)

parser.add_argument(
    "-v",
    action="count",
    default=0,
    dest="verbosity_level",
    help="""
        Use -v to print out information and commands as the test runs and -vv to
        show outputs of commands as the tests run. Note, this information is always
        displayed if there is a test failure""",
)

if __name__ == "__main__":
    sys.exit(main(parser.parse_args()))
