#!/usr/bin/env python3
import argparse
import os
import subprocess


def get_last_line(path):
    last_line = None
    with open(path, "r") as f:
        for line in f:
            if len(line) > 0 and not line.isspace():
                last_line = line
    if last_line is None:
        raise ValueError("the {} file is empty".format(path))
    return last_line.rstrip()


def query_version():
    return get_last_line(
        os.path.join(os.path.dirname(__file__), "..", "VERSION")
    )


def _call(command, fallback_result=None, **kwargs):
    try:
        rslt = subprocess.check_output(command, shell=True, **kwargs)
    except (subprocess.CalledProcessError, OSError):
        return fallback_result
    return rslt.decode().rstrip()  # return as str & remove any trailing '\n'


def query_git(command):
    # historically, we queried whether git exists before executing a command.
    # However, on certain systems this doesn't seem to be adequate. Instead, we
    # just try to execute the command. If it fails, we fall back to "N/A"
    return _call(command, fallback_result="N/A")


choices = {
    "show-version": query_version,
    "git-branch": lambda: query_git("git rev-parse --abbrev-ref HEAD"),
    "git-revision": lambda: query_git("git rev-parse HEAD"),
}

parser = argparse.ArgumentParser("query version information")
parser.add_argument(
    "directive",
    choices=list(choices),
    help="specifies the information to check",
)

if __name__ == "__main__":
    args = parser.parse_args()
    result = choices[args.directive]()
    print(result)
