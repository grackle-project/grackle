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
    

def configure_file(lines, variable_map, out_fname, literal_linenos):
    """
    Writes a new file to out_fname, line-by-line, while performing variable
    substituions
    """

    used_variable_set = set()
    out_f = open(out_fname, 'w')
    err_msg = None

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

    for line_num, line in enumerate(lines):
        # make sure to drop any trailing '\n'
        assert line[-1] == '\n', "sanity check!"
        line = line[:-1]
        if line_num in literal_linenos:
            subbed = line
        else:
            match_count = 0
            subbed = _PATTERN.sub(replace,line)
        out_f.write(subbed)
        out_f.write('\n')
        if err_msg is not None:
            out_f.close()
            os.remove(out_fname)
            raise RuntimeError(err_msg)

    unused_variables = used_variable_set.symmetric_difference(variable_map)

    if len(unused_variables) > 0:
        os.remove(out_fname)
        raise RuntimeError("the following variable(s) were specified, but "
                           "were unused: {!r}".format(unused_variables))

def _parse_variables(dict_to_update, var_val_assignment_str_l,
                     val_kind = 'literal'):
    assert val_kind in ['literal', 'file-path-escaped-contents',
                        'file-path-literal-contents']
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

        if val_kind != 'literal': # val_kind is some kind of file path
            path = value
            if not os.path.isfile(path):
                raise RuntimeError(
                    ("error while trying to associate the contents of the file "
                     "at {!r} with the {!r} variable: no such file exists"
                     ).format(path, var_name))
            with open(value, "r") as f:
                if val_kind == 'file-path-escaped-contents':
                    # we generally treat the characters in the file as literals
                    # -> we do need to make a point of properly escaping the
                    #    newline characters
                    assert os.linesep == '\n' # implicit assumption
                    value = f.read().replace(os.linesep, r'\n')
                else: # val_kind == 'file-path-literal-contents'
                    value = f.read()

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
                     val_kind = 'literal')
    _parse_variables(variable_map, args.variable_use_file_contents,
                     val_kind = 'file-path-escaped-contents')
    _parse_variables(variable_map, args.variable_use_literal_file_contents,
                     val_kind = 'file-path-literal-contents')

    literal_linenos = set()
    if args.literal_linenos is not None:
        literal_linenos = set(args.literal_linenos)
    # use variable_map to actually create the output file
    with open(args.input, 'r') as f_input:
        line_iterator = iter(f_input)
        configure_file(lines = line_iterator,
                       variable_map = variable_map,
                       out_fname = out_fname,
                       literal_linenos=literal_linenos)

    return 0

parser = argparse.ArgumentParser(description='Configure template files.')
parser.add_argument(
    '--variable-use-file-contents',  action = 'append', default = [],
    metavar = 'VAR=path/to/file',
    help = ("associates the (possibly multi-line) contents contained by the "
            "specified file with VAR. This replaces each newline character "
            "with the pair of characters \"\\n\". This is useful if the "
            "contents represent a string to be printed")
)

parser.add_argument(
    '--variable-use-literal-file-contents',  action = 'append', default = [],
    metavar = 'VAR=path/to/file',
    help = ("associates the (possibly multi-line) contents contained by the "
            "specified file with VAR. This does NOT escape newline characters.")
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
parser.add_argument(
    "--literal-linenos", nargs="*", type=int,
    help = "line numbers corresponding to lines that are treated as literals"
)

if __name__ == '__main__':
    sys.exit(main(parser.parse_args()))
