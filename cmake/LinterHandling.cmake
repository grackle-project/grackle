# This module defines functionality to help configure cmake to use sanitizers

# Clang-Tidy Support
# ==================
#
# Refer to the website documentation for a more detailed overview of what
# clang-tidy is and how you use it.
#
# Goals of this functionality
# ---------------------------
# 1. Make it easier for non-CMake experts to use clang-tidy (in a totally
#    "opt-in" manner so that CMake experts can do things in the standard way)
# 2. Backport a change (introduced in CMake 3.25) that alters the arguments
#    passed to CMake when the `-p` flag is detected. It's important to pass
#    this flag to work around a well known bug (described further down below)
#
# How to use this functionality
# -----------------------------
# Execute `convenience_prepare_clang_tidy()` macro at the top-level of the
# CMake project BEFORE declaring any targets (the macro's docstring provides
# further context)


# a helper function that constructs a clang-tidy command line that effectively
# backports behavior introduced CMake 3.25 when the -p flag is passed
#
# Desired Behavior:
# - Ordinarily, CMake invokes clang-tidy by appending some extra arguments
#   related to the source-file, the 2-character `--` argument, and then more
#   general compilation options after the `--` argument
# - starting in CMake 3.25, when the `-p` argument is passed, all arguments
#   including and following the `--` argument are omitted
#
# Arguments:
# - `TIDY_COMMAND_LINE` holds a value that you should be able to directly store
#   in a target's `<LANG>_CLANG_TIDY` property (in CMake 25 or later). The docs
#   describe this as "a semicolon-separated list containing a command line for
#   the clang-tidy tool." (The list's 1st element specifies executed command &
#   subsequent items are extra args).
# - `OUT_VAR` is the variable-name where the adjusted "command-line for
#   clang-tidy is written"
#
# Side-Effects:
# To achieve the desired behavior, we create a wrapper script that scans the
# arguments, makes any adjustments, and forwards them onto clang-tidy
function(_preptidy_backport_cmake325_cli_handling TIDY_COMMAND_LINE OUT_VAR)

  set(wrapper_script ${CMAKE_CURRENT_BINARY_DIR}/clang-tidy-wrapper)

  # at the moment, the script expects that the very first argument is the path
  # to the clang-tidy. We could maybe revisit this in the future once we
  # require CMake >= 3.18 and can use `file(CONFIGURE ...)`
  file(WRITE ${wrapper_script} [=[
#!/usr/bin/env python3
import sys, subprocess

def filter_args(args):
    for i, v in enumerate(args[1:], start=1): # arg[0] is path/to/clang-tidy
        if v == "--":
            break
        elif v == "-p":
            try:
                return args[:args.index("--",i+1)]
            except ValueError:
                break
    return args

unfiltered_args = sys.argv[1:] # sys.argv[0] is equiv to __file__
sys.exit(subprocess.call(filter_args(unfiltered_args)))
]=])

  file(CHMOD ${wrapper_script}
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ)

  set("${OUT_VAR}" "${wrapper_script};${TIDY_COMMAND_LINE}" PARENT_SCOPE)

endfunction()

# this is a helper function used within convenience_prepare_clang_tidy
# -> this function manipulates actually uses GRACKLE_DEV_CLANG_TIDY variable to
#    set CMAKE_CXX_CLANG_TIDY (see the docstring of
#    convenience_prepare_clang_tidy for details)
function(_preptidy_helper)

  if(GRACKLE_IS_TOP_LEVEL AND
      (${CMAKE_MINIMUM_REQUIRED_VERSION} VERSION_GREATER_EQUAL "3.25"))
    message(FATAL_ERROR
      "this is a reminder, now that the minimum required CMake version is at "
      "least 3.25, that the CMake logic can be simplified. Feel free to remove "
      "this logic after updating the issue-tracker")
  endif()

  list(FIND GRACKLE_DEV_CLANG_TIDY "-p" index)
  if (index GREATER_EQUAL 1)
    message(FATAL_ERROR
      "if you know enough to pass the -p flag to GRACKLE_DEV_CLANG_TIDY, then you "
      "should set the CMAKE_CXX_CLANG_TIDY variable, directly"
    )
  endif()

  # determine the absolute path to clang-tidy
  # -> for CMake >= 3.25, this is a convenience to help prevent people from
  #    "shooting themselves in the foot"
  # -> in early versions of CMake, this *may* be necessary so that our
  #    backported functionality works properly (but I haven't checked carefully)
  list(GET GRACKLE_DEV_CLANG_TIDY 0 program_name)

  if (GRACKLE_DEV_CLANG_TIDY STREQUAL "")
    message(FATAL_ERROR "GRACKLE_DEV_CLANG_TIDY holds an empty string")
  elseif (program_name MATCHES "/")
    set(path ${program_name})
    if (NOT IS_ABSOLUTE "${path}")
      message(FATAL_ERROR
        "GRACKLE_DEV_CLANG_TIDY specifies a path, `${path}`, that isn't absolute")
    elseif (NOT EXISTS "${path}")
      message(FATAL_ERROR
        "GRACKLE_DEV_CLANG_TIDY specifies a path, `${path}`, that doesn't exist")
    endif()
  else()
    message("searching for ${program_name} since it doesn't seem to be a path")
    # the output variable is a cache variable
    find_program(_preptidy_CLANG_TIDY_PATH ${program_name})
    if (_preptidy_CLANG_TIDY_PATH) # the search was succesful
      set(path ${_preptidy_CLANG_TIDY_PATH})
    else() # _preptidy_CLANG_TIDY_PATH holds value ending in NOTFOUND
      message(FATAL_ERROR
        "Could not find a program named \"${program_name}\"")
    endif()
  endif()

  # rebuild a version of GRACKLE_DEV_CLANG_TIDY that uses path
  set(coerced_cmd "${path}")
  list(LENGTH GRACKLE_DEV_CLANG_TIDY n_components)
  if (n_components GREATER 1)
    list(SUBLIST GRACKLE_DEV_CLANG_TIDY 1 -1 specified_args)
    set(coerced_cmd "${coerced_cmd};${specified_args}")
  endif()

  # adjust coerced_cmd for older CMake versions
  if(CMAKE_VERSION VERSION_LESS "3.25")
    _preptidy_backport_cmake325_cli_handling("${coerced_cmd}" tmp)
    set(coerced_cmd "${tmp}")
  endif()

  # instruct clang-tidy to use the compilation database that CMake generates
  # -> this avoids lots of issues with clang-tidy mis-identifying toolchain
  #    headings
  set(CMAKE_CXX_CLANG_TIDY "${coerced_cmd};-p;${CMAKE_BINARY_DIR}" PARENT_SCOPE)

endfunction()


# Setup CMAKE_CXX_CLANG_TIDY variable from the GRACKLE_DEV_CLANG_TIDY variable
# (does nothing if the GRACKLE_DEV_CLANG_TIDY variable isn't defined).
# 
# This macro is intended purely as a convenienced to help properly configure
# clang-tidy for contributors (who aren't necessarily CMake experts) to use the
# tool locally. More experienced users are welcomed to directly set the
# CMAKE_CXX_CLANG_TIDY variable.
#
# Specifically this convenience method ensures that CMake invokes clang-tidy
# using the `-p` flag. This flag instructs clang-tidy to use the compilation
# database that CMake generates as part of the build. This avoids a bug where
# clang-tidy can misidentify the toolchain headers being in a given build
#
# How this should be called
# =========================
# This macro should be called at the top level scope, after defining the
# GRACKLE_IS_TOP_LEVEL variable and ensuring that the
# CMAKE_EXPORT_COMPILE_COMMANDS variable is defined
#
# About the GRACKLE_DEV_CLANG_TIDY argument
# =========================================
# (Aside: keep this description synchronized with the docs)
#
# This variable should specify a semicolon-delimited list, where the first
# list element specifies the name of the `clang-tidy` executable (CMake searches
# for the executable in all the standard places) or the absolute path to the
# executable. Any subsequent list elements specify command-line arguments are
# forwarded to the executable (e.g. a common choice is `-warnings-as-errors=*`)
#
# Since this is a developer variable, (and NOT a configuration variable for
# users), I think we are safe to adjust the use/behavior of this variable
macro(convenience_prepare_clang_tidy)
  # do some basic sanity checks to verify that edits to the CMakeLists.txt file
  # did not break any assumptions made by this macro
  if (DEFINED CMAKE_CURRENT_FUNCTION)
    message(fatal_error
      "convenience_prepare_clang_tidy macro: called from within a function")
  elseif(NOT DEFINED GRACKLE_IS_TOP_LEVEL)
    message(fatal_error
      "convenience_prepare_clang_tidy macro: called without defining "
      "GRACKLE_IS_TOP_LEVEL")
  elseif(GRACKLE_IS_TOP_LEVEL AND NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
    message(fatal_error
      "convenience_prepare_clang_tidy macro: called without defining "
      "CMAKE_EXPORT_COMPILE_COMMANDS (and Grackle is the top-level project)")
  endif()

  # onto the main logic:
  if (DEFINED GRACKLE_DEV_CLANG_TIDY)
    # check if the end-user made any mistakes
    if ((DEFINED CMAKE_C_CLANG_TIDY) OR (DEFINED CMAKE_CXX_CLANG_TIDY))
      message(FATAL_ERROR
        "It's an error to define GRACKLE_DEV_CLANG_TIDY while also defining "
        "CMAKE_C_CLANG_TIDY or CMAKE_CXX_CLANG_TIDY")
    elseif (NOT "${GRACKLE_IS_TOP_LEVEL}")
      message(FATAL_ERROR
        "It's an error to define GRACKLE_DEV_CLANG_TIDY when Grackle isn't the "
        "top-level project")
    elseif (NOT "${CMAKE_EXPORT_COMPILE_COMMANDS}")
      message(FATAL_ERROR
        "It's an error to define GRACKLE_DEV_CLANG_TIDY when "
        "CMAKE_EXPORT_COMPILE_COMMANDS is turned off")
    elseif (NOT DEFINED CACHE{GRACKLE_DEV_CLANG_TIDY})
      # it's important to enforce that GRACKLE_DEV_CLANG_TIDY is a cache entry
      # so users get consistent behavior if cmake runs again between rebuilds
      message(FATAL_ERROR
        "somehow the GRACKLE_DEV_CLANG_TIDY variable isn't a cache variable. "
        "This shouldn't be possible if you used set the value on the "
        "command-line with -DGRACKLE_DEV_CLANG_TIDY=<value>"
      )
    else()
      _preptidy_helper()
    endif()
  endif()

endmacro()

