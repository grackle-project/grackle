#!/usr/bin/env python3

########################################################################
#
# standalone python script that is used to help check code examples
# - this is NOT a direct part of the pytest framework
# - for compatability, we explicitly choose to:
#   - only use modules from the standard library
#   - only use features from python 3.7 or newer
#
# Copyright (c) 2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import argparse
import functools
import json
import math
import os
import re
import stat
import subprocess
import sys
import textwrap

# Common machinery used for multiple subcommands
# ----------------------------------------------

rfields = (
    "cooling_time",
    "dust_temperature",
    "gamma",
    "pressure",
    "temperature"
)

_LINE_REGEX = re.compile(rf"^ ?({'|'.join(rfields)}) = (.*)$")

# following line comes from the re module's docs for simulating scanf:
_FLT_PATTERN =r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
_FIELD_VAL_REGEX = re.compile(
    rf"^\s*(?P<value>{_FLT_PATTERN})(\s*|\s+(?P<units>.+))$"
)

def parse_output(ostr):
    """Parse the output of a code-example"""
    results = {field: None for field in rfields}

    ostr = ostr.decode("utf8") if isinstance(ostr, bytes) else ostr
    lines = ostr.split("\n")
    for line in lines:
        match = _LINE_REGEX.match(line)
        if match is None:
            continue
        field = match.group(1)
        rside = match.group(2)

        match = _FIELD_VAL_REGEX.match(rside)
        if match is None:
            raise RuntimeError(
                f"Cannot grab field values from line: {line}.")

        d = match.groupdict()
        val = float(d['value'])
        units = 'dimensionless' if d['units'] is None else d['units']

        if results[field] is not None:
            raise RuntimeError(
                "code example output multiple lines specifying the value for "
                f"{field}.")
        results[field] = (val, units)

    return results

def run_and_parse(full_path, timeout, env=None, exec_dir=None,
                  announce_execution=True):
    """Executes a code-example and returns the parsed output"""
    timeout = None if timeout==0 else timeout

    if exec_dir is None:
        cwd = os.path.dirname(full_path)
        cwd = None if cwd == '.' else cwd
        exec_path = os.path.join('.', os.path.basename(full_path))
    else:
        exec_path = os.path.abspath(full_path)
        cwd = exec_dir

    if announce_execution:
        print(
            f" excuted cmd:  {exec_path}",
            f" working dir:  {cwd}",
            f" timeout:      {timeout}",
            f" env-override: {env}",
            sep="\n"
        )
    proc = subprocess.run(
        exec_path, timeout=timeout, cwd=cwd, env=env, capture_output=True
    )
    if proc.returncode == 0:
        log = proc.stdout
    else:
        raise RuntimeError(
            f"Command {command} failed with return value "
            f"{proc.returncode} and the following stderr output "
            f"{proc.stderr}")

    return parse_output(log)

# generic parsering logic
# -----------------------

_exec_stat_mask = stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH
def _is_executable(path):
    return (os.stat(path).st_mode & _exec_stat_mask) > 0

def _path_to_existing_file(s, req_executable=False):
    # intended to be used as a type kwarg of add_argument
    if not os.path.isfile(s):
        raise ValueError(f"no file called `{s}`")
    elif req_executable and not _is_executable(s):
        raise ValueError(f"`{s}` is not executable")
    return s

def _path_to_existing_dir(s):
    if not os.path.isdir(s):
        raise ValueError(f"no directory called `{s}`")
    return s

def _nonnegative_finite_float(s):
    # intended to be used as a type kwarg of add_argument
    val = float(s)
    if (val < 0) or not math.isfinite(val):
        raise ValueError(f'`{s}` is not a non-negative finite float')
    return val

def _add_exec_timeout_flag(parser):
    parser.add_argument(
        "--exec-timeout",
        metavar="DURATION",
        default=0.0,
        required=False,
        type=_nonnegative_finite_float,
        help=(
            "when execution of a code-example exceeds the specified DURATION "
            "(a floating point number with units of seconds), the program "
            "exits with an error code denoting failure"
        )
    )

def _add_exec_dir_flag(parser, required=False):
    extra_help = ""
    if not required:
        extra_help = (
            " When not specified, the default behavior is to execute from the "
            "directory where the binary is defined"
        )

    parser.add_argument(
        "--exec-dir",
        required=required,
        type=_path_to_existing_dir,
        help=(
            "directory where the code-examples binary are executed from."
            + extra_help
        )
    )

# genjson subcommand
# ------------------

def _add_genjson_subcommand(subparsers):
    genjson_parser = subparsers.add_parser(
        "genjson",
        help="generate a json file that stores the results of a code-example"
    )
    genjson_parser.add_argument(
        "--target",
        required=True,
        type=functools.partial(_path_to_existing_file, req_executable=True),
        help=(
            "path to the compiled executable that will be executed. Relative "
            "path are relative to caller's directory (i.e. it is independent "
            "of --exec-dir)"
        )
    )
    genjson_parser.add_argument(
        "--output", required=True, help="path to the output json file"
    )
    _add_exec_timeout_flag(genjson_parser)
    _add_exec_dir_flag(genjson_parser, required=False)
    genjson_parser.add_argument('-v', '--verbose', action='store_true')
    genjson_parser.set_defaults(func = main_genjson)


def main_genjson(args):
    """Implements the `genjson` subcommand."""
    print(f'executing `{args.target}`')
    rslt = run_and_parse(
        args.target, timeout=args.exec_timeout, exec_dir=args.exec_dir
    )
    if args.verbose:
        print(rslt)

    print(f"recording results to `{args.output}`")
    with open(args.output, mode="w") as f:
        json.dump(obj=rslt, fp=f, indent=4)

# cmp subcommand
# --------------

def _add_cmp_subcommand(subparsers):
    cmp_parser = subparsers.add_parser(
        "cmp", help="compare the results of 2 code-example executables"
    )
    for flag,descr in [("--ref", "used as a reference"),
                       ("--target", "tested")]:
        cmp_parser.add_argument(
            flag,
            required=True,
            type=_path_to_existing_file,
            help=(
                "path to the compiled code-executable to be {descr} (OR to a "
                "json file recorded the results of a prior execution). "
                "Relative path are relative to caller's directory (i.e. it is "
                "independent of --exec-dir)"
            )
        )

    for flag, kind in [("--atol","absolute"), ("--rtol","relative")]:
        cmp_parser.add_argument(
            flag,
            default=0.0,
            type=float,
            required=False,
            help=f"specifies the {kind} tolerance used for the comparison"
        )
    _add_exec_timeout_flag(cmp_parser)
    _add_exec_dir_flag(cmp_parser, required=False)
    cmp_parser.add_argument('-v', '--verbose', action='store_true')
    cmp_parser.set_defaults(func = main_cmp)

def main_cmp(args):
    """Implements the `cmp` subcommand."""
    result_l = []

    for attr in ['ref', 'target']:
        path = getattr(args, attr)
        if _is_executable(path):
            print(f'{attr}: executing {path}')
            rslt = run_and_parse(
                path, timeout=args.exec_timeout, exec_dir=args.exec_dir
            )
        else:
            print(f'{attr}: reading {path}')
            with open(path, 'r') as f:
                rslt = {k: tuple(v) for k,v in json.load(f).items() }

        result_l.append(rslt)

    ref_rslt, target_rslt = result_l
    atol, rtol = args.atol, args.rtol
    print(f"Comparing Results: atol = {atol}, rtol = {rtol}")

    successful_check = isclose_resultpacks(
        actual=target_rslt,
        reference=ref_rslt,
        atol=args.atol,
        rtol=args.rtol,
        verbose=args.verbose
    )
    if successful_check:
        print("SUCCESS")
        return 0
    return 1

class CmpErrorRecorder:
    # as we perform a comparison, instances of this class record errors
    def __init__(self): self.log = []
    def record_err(self, s): self.log.append(s)
    def any_err(self): return len(self.log) > 0
    def __iter__(self): return iter(self.log)

def _consistent_keys(actual, reference, err_recorder):
    # returns iterable for the keys shared by actual and reference
    refkeys = reference.keys()
    refkey_set = set(refkeys)
    mismatch_keys = refkey_set.symmetric_difference(actual.keys())
    if len(mismatch_keys) == 0:
        return actual.keys()
    else:
        shared_keys = list(refkey_set.intersection(actual.keys()))
        extra_ref, extra_actual = [], []
        for k in mismatch_keys:
            if k in refkeys:
                extra_ref.append(k)
            else:
                extra_actual.append(k)
        err_recorder.record_err(textwrap.dedent(f"""\
            Key Mismatch Error: the following are extra (unmatched) keys
              actual:    {extra_actual!r}
              reference: {extra_ref!r}"""
        ))
        return shared_keys


def isclose_resultpacks(actual, reference, rtol, atol, verbose=False):
    """
    Returns True if entries of the result packs are equal within a tolerance.

    Parameters
    ----------
    actual, reference : mapping
        A mapping holding results obtained in calculations. The keys are the
        names of quantities. Each value is a tuple where the first element is
        a number and the second value is a string (corresponding to a uniy)
    rtol : float
        The relative tolerance parameter
    atol : float
        The absolute tolerance parameter
    verbose : bool
        Print information
    announce_mismatch : bool
        Print detailed error-reports on failure
    """

    if verbose:
        print(f"Comparing Results with atol = {atol}, rtol = {rtol}")

    match_count = 0
    checked_cases = 0
    err_recorder = CmpErrorRecorder()
    keys = _consistent_keys(actual, reference, err_recorder)
    for key in keys:
        checked_cases += 1
        actual_pair, ref_pair = actual[key], reference[key]

        if verbose:
            print(f"  comparing {key}: actual={actual!r}, ref={ref!r}")

        if (ref_pair is None) or (actual_pair is None):
            if actual_pair != ref_pair:
                err_recorder.record_err(textwrap.dedent(f"""\
                    {key!r}: one of the values is None
                      actual: {actual_pair!r}
                      ref: {ref_pair!r}"""
                ))
            else:
                match_count += 1
            continue

        # we use numpy's convention of asymmetric relative differences
        ref_val, ref_units = ref_pair
        actual_val, actual_units = actual_pair
        diff = actual_val - ref_val

        if ref_units != actual_units:
            err_recorder.record_err(textwrap.dedent(f"""\
                {key!r}: Unit mismatch!
                  units: actual = {actual_units!r}, ref = {ref_units!r}
                  values: actual = {actual_val!r}, ref = {ref_val!r}"""
            ))
        elif (abs(diff) > (atol + rtol * abs(ref_val))):
            rdiff_str = f'{diff/ref_val:g}' if (ref_val != 0) else 'undefined'
            err_recorder.record_err(textwrap.dedent(f"""\
                {key!r}: not equal to specified tolerance
                  relative Difference: {rdiff_str}
                  absolute Difference: {abs(diff):g}
                  values: actual = {actual_val!r}, ref = {ref_val!r}
                  units: {ref_units!r}"""
            ))
        else:
            match_count+=1

    if not err_recorder.any_err():
        return True

    msg = textwrap.dedent(f"""
        FAIL: Result-packs not equal to tolerance rtol={rtol}, atol={atol}
          Context:
            compared results for the following keys:
              {repr(list(keys))[1:-1]}
            {match_count} of these {checked_cases} results matched"""
    )
    print(
        msg, *(textwrap.indent(e,prefix='  ') for e in err_recorder), sep='\n'
    )
    return False


# datacheck subcommand
# --------------------

def _add_datacheck_subcommand(subparsers):
    datacheck_parser = subparsers.add_parser(
        "datacheck",
        help=(
            "checks whether code examples' hardcoded assumption about the path "
            "to the data file is satisfied"
        )
    )
    datacheck_parser.add_argument(
        "--assumed-data-path",
        required=True,
        help=(
            "Hardcoded path to the data file in the code-examples (relative "
            "to the location in --exec-dir"
        )
    )
    _add_exec_dir_flag(datacheck_parser, required=True)
    datacheck_parser.set_defaults(func = main_datacheck)

def main_datacheck(args):
    """Implements the `datacheck` subcommand."""
    if os.path.isabs(args.assumed_data_path):
        path = args.assumed_data_path
    else:
        path = os.path.join(args.exec_dir, args.assumed_data_path)
    if os.path.isfile(path):
        print("File exists!")
        return 0
    else:
        print(f"No file at {path}")
        return 1

# Tie everything together
# -----------------------

def main(args):
    return args.func(args)

parser = argparse.ArgumentParser(
    prog="code_example_checker",
    description="tool that helps check correctness of code examples",
)

subparsers = parser.add_subparsers(help="subcommand help", required=True)
_add_genjson_subcommand(subparsers)
_add_cmp_subcommand(subparsers)
_add_datacheck_subcommand(subparsers)

if __name__ == '__main__':
    sys.exit(main(parser.parse_args()))
