#!/usr/bin/env python3

# for the sake of portability, we try to avoid fstrings

import argparse
import os
import re
import string
import sys

_MAX_VARNAME_SIZE = 256

def is_valid_varname(s, start = None, stop = None):
    return re.fullmatch("\w+", s[slice(start, stop)]) is not None

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

    length = len(line)
    pos = 0

    while pos < length:
        # Step 1: Scan until @ symbol and copy all characters from `pos` up to
        #         `leading` to out_f
        leading = line.find('@', pos)
        if leading == pos:
            pass # do nothing - no characters to write
        elif leading == -1:
            # there aren’t any variable substitutions in line[pos:]
            out_f.write(line[pos:])
            break
        else:
            out_f.write(line[pos:leading])

        # The remainder of this loop handles variable substitution

        # Step 2: scan until trailing @ symbol
        trailing = line.find('@', leading + 1,
                             _MAX_VARNAME_SIZE + 2 + leading)
        if trailing == -1:
            return ("the '@' in column {} of line {} is isolated, the "
                    "specified variable isn’t properly terminated, OR the "
                    "variable name is too long").format(leading, line_num)
        elif trailing == (leading+1):
            return "\"@@\" appears at column {} of line {}".format(leading,
                                                                   line_num)
        elif not is_valid_varname(line, leading + 1, trailing):
            return ("the string enclosed by the '@' characters, starting at "
                    "column {} of line {} doesn’t specify an allowed variable "
                    "name").format(leading, line_num)

        # Step 3: write the variable value to the output file
        varname = line[leading+1:trailing]
        if varname not in variable_map:
            return ("the variable {} (specified by a string enclosed by the "
                    "'@' characters at column {} of line {}) doesn't have an "
                    "associated value").format(varname, leading, line_num)
        used_variable_set.add(varname)
        out_f.write(variable_map[varname])

        # Step 4: update pos:
        pos = trailing+1

    # after leaving the loop, write a \n
    out_f.write('\n')
    return None


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
