#!/bin/sh

# This script intentionally tries to just use features offered in the posix
# shell to maximize compatability (i.e. bashisms are avoided)
# -> The primary exception is the local keyword, which is not part of POSIX,
#    but is supported by nearly all POSIX-compatible shells
# -> If we aren't worried about compatability, we might be better off rewriting
#    this in python (could support 2 or 3)
#
# Most of the code in this file was ripped out of the build-system


prog="$0"
usage_and_exit(){
    echo ""
    echo "Utility program designed to help with building Grackle"
    echo ""
    echo "Usage:"
    printf "    %s show-version -version-file <path> \n" "$prog"
    printf "    %s generate-funcs -version-file <path> \n" "$prog"
    printf "        -config-file <path> -flag-file <path> -out <path>\n"
    printf "    %s generate-header -precision (32 | 64) -out <path>\n" "$prog"
    echo ""
    echo "Flags:"
    echo "    -version-file <path>"
    echo "        points to the file containing the version number"
    echo "    -config-file <path>"
    echo "        points to the file containing the config information."
    echo "    -flag-file <path>"
    echo "        points to the file holding compilation flags."
    echo "    -out <path>"
    echo "        specifies the path to the file where output is written"
    echo "    -precision (32 | 64)"
    echo "        specifies the floating point precision. Either 32 or 64."
    if [ $# -gt 0 ]; then
        echo ""
        echo "$@"
    fi
    exit 1
}

# useful functions:

in_list(){ # check whether $1 is contained within $2, a ;-separated list
    local _item="$1"
    local _list="$2"

    local _elem=
    while [ "${_list}" != "${_elem}" ] ;do
        _elem=${_list%%;*} # extract substring up to (but not including) ';'
        _list="${_list#$_elem;}" # set _list to everything after 1st occurance
                                 # of _elem + ';'

        if [ "${_item}" = "${_elem}" ]; then echo "true" && return 0; fi
    done
    echo "false"
}

assert_opt_has_mode(){ # abort with error message if the mode global variable
                       # is not contained within the ;-separated list
    local _opt="$1"
    local _list="$2"
    if [ $(in_list "${mode}" "${_list}") != "true" ]; then 
        usage_and_exit "ERR: \"${mode}\" mode incompatible with \"${_opt}\""
    fi
}

is_help_flag(){ # prints true if $1 is recognized as a help flag
    case "${1}" in
        "-?"|"-h"|"-help"|"--help") echo "true"; return 0;;
        *) echo "false"; return 0;;
    esac
}

# ==================================================================
# PARSE Command line args
# ==================================================================

# parse the mode argument
if [ $# -eq 0 ]; then
    usage_and_exit "ERR: no arguments"
elif [ $(is_help_flag "$1") = "true" ]; then
    usage_and_exit
elif [ $(in_list "$1" "show-version;generate-funcs;generate-header") = "true" ]
then
    mode="$1"
    shift
else
    usage_and_exit "ERR: $0 executed with unexpected mode: $1"
fi

# parse the remaining arguments
configfile=
flagfile=
versionfile=
outfile=
precision=

while [ $# -gt 0 ]
do
    cur_flag="${1}"
    shift
    case "${cur_flag}" in

        "-version-file")
            assert_opt_has_mode "${cur_flag}" "show-version;generate-funcs"
            versionfile="${1}"
            shift
            ;;

        "-config-file")
            assert_opt_has_mode "${cur_flag}" "generate-funcs"
            configfile="${1}"
            shift
            ;;

        "-flag-file")
            assert_opt_has_mode "${cur_flag}" "generate-funcs"
            flagfile="${1}"
            shift
            ;;

        "-out")
            assert_opt_has_mode "${cur_flag}" "generate-funcs;generate-header"
            outfile="${1}"
            shift
            ;;

        "-precision")
            assert_opt_has_mode "${cur_flag}" "generate-header"
            precision="${1}"
            shift
            ;;

        *)
            if [ $(is_help_flag "${cur_flag}") = "true" ]; then
                usage_and_exit
            else
                usage_and_exit "ERR: Unrecognized parameter: ${cur_flag}"
            fi
            ;;
    esac
done

# ==================================================================
# DEFINE functions that do heavy-lifting in each mode
# ==================================================================

VERSION_NUM=
GIT_BRANCH=
GIT_REVISION=
fetch_version_info(){
    if [ -z "$versionfile" ]; then
        usage_and_exit "ERR: -version-file MUST be specified\n"
    else
        VERSION_NUM="$(tail -1 $versionfile)"
    fi

    # note: `&>` is not portable, so we explicitly redirect stderr
    local USEGIT=$(command -v git > /dev/null 2>&1 &&
                   git status > /dev/null 2>&1 &&
                   echo "1")

    if [ "${USEGIT}" = "1" ]; then
        GIT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
        GIT_REVISION="$(git rev-parse HEAD)"
    else
        GIT_BRANCH="N/A"
        GIT_REVISION="N/A"
    fi
}

show_version(){
    fetch_version_info

    echo "The Grackle Version ${VERSION_NUM}"
    echo "Git Branch   ${GIT_BRANCH}"
    echo "Git Revision ${GIT_REVISION}"
    exit 0
}

generate_funcs(){

    # confirm we have required parameters:
    if [ -z "$versionfile" ]; then
        usage_and_exit "ERR: -version-file MUST be specified\n"
    elif [ -z "$configfile" ]; then
        usage_and_exit "ERR: -config-file MUST be specified\n"
    elif [ -z "$flagfile" ]; then
        usage_and_exit "ERR: -flag-file MUST be specified\n"
    elif [ -z "$outfile" ]; then
        usage_and_exit "ERR: -out-file MUST be specified\n"
    fi

    fetch_version_info # fetch relevant version info

    # delete outfile if it already exists
    if [ -f $outfile ]; then rm $outfile; fi

    # generate the output file:
    echo '// this file was automatically generated'                 >> $outfile
    echo ''                                                         >> $outfile
    echo '#include <stdio.h>'                                       >> $outfile
    echo '#include "grackle_types.h"'                               >> $outfile
    echo ''                                                         >> $outfile
    echo 'grackle_version get_grackle_version(void) {'              >> $outfile
    echo '  grackle_version out;'                                   >> $outfile
    echo "  out.version = \"${VERSION_NUM}\";"                      >> $outfile
    echo "  out.branch = \"${GIT_BRANCH}\";"                        >> $outfile
    echo "  out.revision = \"${GIT_REVISION}\";"                    >> $outfile
    echo '  return out;'                                            >> $outfile
    echo '}'                                                        >> $outfile
    echo ''                                                         >> $outfile
    echo 'void auto_show_flags(FILE *fp) {'                         >> $outfile
    awk '{print "  fprintf (fp,\""$$0"\\n\");"}'        < $flagfile >> $outfile
    echo "}"                                                        >> $outfile
    echo ''                                                         >> $outfile
    echo 'void auto_show_config(FILE *fp) {'                        >> $outfile
    awk '{print "  fprintf (fp,\""$$0"\\n\");"}'      < $configfile >> $outfile
    echo "}"                                                        >> $outfile

    exit 0
}

generate_header(){
    # confirm we have required parameters & determine macro_name:
    local macro_name=""

    if [ -z "$outfile" ]; then
        usage_and_exit "ERR: -out-file MUST be specified\n"
    elif [ -z "${precision}" ]; then
        usage_and_exit "ERR: -precision MUST be specified\n"
    elif [ "${precision}" = "32" ]; then
        macro_name="GRACKLE_FLOAT_4"
    elif [ "${precision}" = "64" ]; then
        macro_name="GRACKLE_FLOAT_8"
    else
        usage_and_exit "ERR: \"-precision\" must be assigned 32 or 64"
    fi
    
    # delete outfile if it already exists
    if [ -f $outfile ]; then rm $outfile; fi

    # generate the output file:

    # skip the following line since the header may be included in fortran files
    #echo '// this file was automatically generated'                >> $outfile

    echo '#ifndef __GRACKLE_FLOAT_H__'                              >> $outfile
    echo '#define __GRACKLE_FLOAT_H__'                              >> $outfile
    echo "#define ${macro_name}"                                    >> $outfile
    echo '#endif'                                                   >> $outfile

    exit 0
}

# ==================================
# DISPATCH to the relevant functions
# ==================================

if [ "$mode" = "show-version" ]; then
    show_version
elif [ "$mode" = "generate-funcs" ]; then
    generate_funcs
elif [ "$mode" = "generate-header" ]; then
    generate_header
else
    usage_and_exit "ERR: Unknown mode: ${mode}. How did we get here?"
fi
