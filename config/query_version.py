#!/usr/bin/env python3
import argparse, os, subprocess

def get_last_line(path):
    last_line = None
    with open(path, 'r') as f:
        for line in filter(lambda l: len(l) > 0 and not l.isspace(), f):
            last_line = line
    if last_line is None:
        raise ValueError("the {} file is empty".format(path))
    return last_line.rstrip()

def query_version():
    return get_last_line(os.path.join(os.path.dirname(__file__), '../VERSION'))

def _call(command, **kwargs):
    rslt = subprocess.check_output(command, shell = True, **kwargs)
    return rslt.decode().rstrip()  # return as str & remove any trailing '\n'

def query_git(command):
    # note: we explicitly redirect stderr since `&>` is not portable 
    git_is_installed = _call('command -v git > /dev/null 2>&1 && '
                             'git status > /dev/null 2>&1 && '
                             'echo "1"') == "1"
    return _call(command) if git_is_installed else "N/A"

choices = {"show-version" : query_version,
           "git-branch" : lambda: query_git("git rev-parse --abbrev-ref HEAD"),
           "git-revision" : lambda: query_git("git rev-parse HEAD")}

parser = argparse.ArgumentParser("query version information")
parser.add_argument('directive', choices = list(choices),
                    help = "specifies the information to check")

if __name__ == '__main__':
    args = parser.parse_args()
    result = choices[args.directive]()
    print(result)
