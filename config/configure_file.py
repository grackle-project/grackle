#!/usr/bin/env python3

# for the sake of portability, we try to avoid fstrings

import argparse
import os
import re
import string
import sys

_MAX_VARNAME_SIZE = 256
_VALID_VARNAME_STR = '\\w{{1,{}}}'.format(_MAX_VARNAME_SIZE)
_PATTERN = re.compile(r'(@{}@)|(@[^\s@]*@?)'.format(_VALID_VARNAME_STR))
_ERR_MSG_TEMPLATE = (
    "{!r}, the string starting with occurence number {} of the '@' character "
    "on line number {} doesn't specify a valid variable name. "
    "A valid variable name is string enclosed by 2 '@' symbolds, where the "
    "string and is composed of 1 to {} alphanumeric ASCII characters. An "
    "alphanumeric character is an uppercase or lowercase letter (A-Z or a-z), "
    "a digit (0-9) or an underscore (_)")

def is_valid_varname(s, start = None, stop = None):
    return re.fullmatch(_VALID_VARNAME_STR, s[slice(start, stop)]) is not None

def process_line(line_num, line, out_f, variable_map, used_variable_set):
    """
    Copies the contents of `line` to `out_f` and perform any variable
    substitutions

    Parameters
    ----------
    line_num : int
        Specifies the line number of the current line (used for formatting
        error messages)
    line : str
        A string representing a line from the input file. This should NOT be
        terminated by a '\n'
    out_f
        Represents the output file
    variable_map : dict
        A dict maping names of variables with the associated values
    used_variable_set : set
        A set updated by this function to record the names of all variables for
        which substitutions have been performed.
    Returns
    -------
    out
        If the function is successful, this is `None`. Otherwise, this is a 
        formatted string describing the error message.
    """

    match_count, err_msg = 0, None

    def replace(matchobj):
        nonlocal err_msg, match_count, used_variable_set, variable_map
        if matchobj.lastindex == 1:
            varname = matchobj[1][1:-1]
            if varname in variable_map:
                match_count += 1
                used_variable_set.add(varname)
                return variable_map[varname]
            err_msg = ("the variable {} (specified by a string enclosed by a "
                       "pair of '@' characters on line {}) doesn't have an "
                       "associated value").format(varname, line_num)
        elif err_msg is None:
            err_msg = _ERR_MSG_TEMPLATE.format(
                matchobj[0], 2*match_count+1, line_num, _MAX_VARNAME_SIZE)
        return '-' # denotes bad case

    out_f.write(_PATTERN.sub(replace,line))
    out_f.write('\n')
    return err_msg

def configure_file(lines, variable_map, out_fname):
    """
    Writes a new file to out_fname, line-by-line, while performing variable
    substituions
    """

    used_variable_set = set()
    out_f = open(out_fname, 'w')

    for line_num, line in enumerate(lines):
        # make sure to drop any trailing '\n'
        assert line[-1] == '\n', "sanity check!"
        line = line[:-1]
        rslt = process_line(line_num = line_num, line = line,
                            out_f = out_f, variable_map = variable_map, 
                            used_variable_set = used_variable_set)
        if rslt is not None:
            out_f.close()
            os.remove(out_fname)
            raise RuntimeError(rslt)

    unused_variables = used_variable_set.symmetric_difference(variable_map)

    if len(unused_variables) > 0:
        os.remove(out_fname)
        raise RuntimeError("the following variable(s) were specified, but "
                           "were unused: {!r}".format(unused_variables))

def _parse_variables(dict_to_update, var_val_assignment_str_l,
                     val_is_file_path = False):
    for var_val_assignment_str in var_val_assignment_str_l:
        stripped_str = var_val_assignment_str.strip() # for safety

        # so the the contents should look like "<VAR>=<VAL>"
        # - For now, we don't tolerate any whitespace.
        # - If we revisit this choice:
        #   - we should be mindful of what it would take to actually escape
        #     whitespace on the command line.
        #   - Doing so often involves quotation marks, that are consumed by the
        #     shell-parsing. It might not be intuitive how such quotation marks
        #     affect the output of this script.
        #   - Consquently it may be more trouble than it's worth to support
        #     whitespace

        for character in (string.whitespace + '@'):
            if character in stripped_str:
                raise RuntimeError(
                    ("a variable-value pair, must not contain the '{}' "
                     "character. The character is present in {!r}").format(
                         character, stripped_str))
        if stripped_str.count('=') != 1:
            raise RuntimeError(
                "each variable-value pair, must contain exactly 1 '=' "
                "charater. This isn't true for {!r}".format(stripped_str))
        var_name, value = stripped_str.split('=')

        if not is_valid_varname(var_name):
            raise RuntimeError(
                "{!r} is not a valid variable name".format(var_name))
        elif var_name in dict_to_update:
            raise RuntimeError(
                "the {!r} variable is defined more than once".format(var_name))

        if val_is_file_path:
            path = value
            if not os.path.isfile(path):
                raise RuntimeError(
                    ("error while trying to associate the contents of the file "
                     "at {!r} with the {!r} variable: no such file exists"
                     ).format(path, var_name))
            with open(value, "r") as f:
                # we generally treat the characters in the file as literals
                # -> we do need to make a point of properly escaping the
                #    newline characters
                assert os.linesep == '\n' # implicit assumption
                value = f.read().replace(os.linesep, r'\n')
        dict_to_update[var_name] = value

def main(args):
    # handle clobber-related logic
    clobber, out_fname = args.clobber, args.output
    if (os.path.isfile(out_fname) and not clobber):
        raise RuntimeError(
            ("A file already exists at {!r}. To overwrite use the --clobber "
             "flag").format(out_fname))

    # fill variable_map with the specified variables and values
    variable_map = {}
    _parse_variables(variable_map, args.variables,
                     val_is_file_path = False)
    _parse_variables(variable_map, args.variable_use_file_contents,
                     val_is_file_path = True)

    # use variable_map to actually create the output file
    with open(args.input, 'r') as f_input:
        line_iterator = iter(f_input)
        configure_file(lines = line_iterator,
                       variable_map = variable_map,
                       out_fname = out_fname)

    return 0

parser = argparse.ArgumentParser(description='Configure template files.')
parser.add_argument(
    '--variable-use-file-contents',  action = 'append', default = [],
    metavar = 'VAR=path/to/file',
    help = ("associates the (possibly multi-line) contents contained by the "
            "specified file with VAR")
)
parser.add_argument(
    "variables", nargs = '*', action = 'store', default = [],
    metavar = 'VAR=VAL',
    help = ("associates the value, VAL, with the specified variable, VAR")
)
parser.add_argument(
    "-i", "--input", required = True, help = "path to input template file"
)
parser.add_argument(
    "-o", "--output", required = True, help = "path to output template file"
)
parser.add_argument(
    "--clobber", action = "store_true",
    help = "overwrite the output file if it already exists"
)

if __name__ == '__main__':
    sys.exit(main(parser.parse_args()))
