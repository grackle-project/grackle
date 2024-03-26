#!/usr/bin/env python3
import argparse, os, subprocess

def _call(command, **kwargs):
    rslt = subprocess.check_output(command, shell = True, **kwargs)
    return rslt.decode()[:-1]  # return as str and remove trailing '\n'

def show_version():
    version_file = os.path.join(os.path.dirname(__file__), '../VERSION')
    print( _call("tail -1 " + version_file) )

def query_git(command):
    # note: we explicitly redirect stderr since `&>` is not portable 
    git_is_installed = _call('command -v git > /dev/null 2>&1 && '
                             'git status > /dev/null 2>&1 && '
                             'echo "1"') == "1"
    result = _call(command) if git_is_installed else "N/A"
    print(result)

choices = {"show-version" : show_version,
           "git-branch" : lambda: query_git("git rev-parse --abbrev-ref HEAD"),
           "git-revision" : lambda: query_git("git rev-parse HEAD")}

parser = argparse.ArgumentParser("query version information")
parser.add_argument('directive', choices = list(choices),
                    help = "specifies the information to check")

if __name__ == '__main__':
    args = parser.parse_args()
    choices[args.directive]()
