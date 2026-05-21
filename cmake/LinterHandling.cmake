# This module defines functionality to help configure cmake to use linters
#
# Current Strategy
# ================
# Our current strategy for supporting linters assumes that developers will
# (directly/indirectly) set the standard CMAKE_<LANG>_<STATIC-ANALYSIS>
# variables (e.g. CMAKE_CXX_CLANG_TIDY or CMAKE_CXX_INCLUDE_WHAT_YOU_USE). For
# short, let's call these the "standard linter variables"
#
# For some context:
# - "standard linter variables" are used to specify the path to binary for the
#   static-analysis tool that the variable is named after and (optionally) some
#   extra args that will be forwarded to the tool
# - we actually provide a convenience function for setting CMAKE_CXX_CLANG_TIDY
#   (see convenience_prepare_clang_tidy)
# - when build-targets are declared, CMake will use the values stored in these
#   variables to set default values of similarly named target-properties (e.g.
#   CXX_CLANG_TIDY or CXX_INCLUDE_WHAT_YOU_USE)
# - During compilation CMake uses these target properties to actually execute
#   the static analysis tools
#
# LinterVar_TmpOverride and LinterVar_Restore
# ===========================================
# There is an unfortunate cross-reaction with using "standard linter variables"
# alongside FetchContent (or git-submodules) to fetch and build dependencies as
# part of the Grackle Build
#
# Explaining the Issue
# --------------------
# Strategies involving FetchContent (or git-submodules) essentially involve
# building dependencies as part of an embedded build. In other words, there is
# almost always a call to `add_subdirectory(<dependency>)` (it happens
# implicitly within `FetchContent_MakeAvailable` and explicitly for
# git-submodules)
#
# Now consider what happens if the are defined before this explicit/implicit
# call to `add_subdirectory`
# - essentially, CMake will try to apply all of the linters to the dependencies
# - this could be problematic because dependencies can definitely fail linter
#   checks
#
# The Solution
# ------------
# To work around this issue, we have created the `LinterVar_TmpOverride` and
# `LinterVar_Restore` variables.
#
# The former command should get called before any calls to
# `FetchContent_MakeAvailable` and the latter should be invoked after all calls
# to `FetchContent_MakeAvailable`. These do something quite crude:
# - `LinterVar_TmpOverride` overrides all "standard linter variables" that are
#   defined and don't hold an empty string. After invoking this comand, any
#   variable expansion ${<standard-linter-variable>} produces an empty string
# - `LinterVar_Restore` reverts all changes made by `LinterVar_TmpOverride`
#
# Alternatives
# ------------
# AFAICT, there isn't really any good robust alternative to this approach.
# - I can think of a few ideas, but they all prevent the direct usage of all
#   built-in "standard-linter-variables" (and we would need to write separate
#   logic to enable each kind of linter)
# - I can think of 2 simplified variants of this approach
#   - starting in CMake 3.25+, we **MIGHT** be able to put the
#     `FetchContent_MakeAvailable` calls within a block() like so:
#     ```
#     block(SCOPE_FOR VARIABLES)
#       set(linter_l "CPPLINT;CLANG_TIDY;CPPCHECK;ICSTAT;INCLUDE_WHAT_YOU_USE")
#       foreach(LINTER IN LISTS linter_list)
#         foreach(LANG IN ITEMS "C" "CXX")
#           set(CMAKE_${LANG}_${LINTER} "")
#         endforeach()
#        endforeach()
#
#        # <inject all of the FetchContent_MakeAvailable logic...>
#
#      endblock()
#      ```
#      We could actually accomplish the exact same effect right now by moving
#      the body of the `block()` into a function and immediately calling the
#      function. While I think this probably works, there's a chance it won't
#      work if CMake assumes that it can propatate info back to the project's
#      primary scope
#    - starting CMake 4.2, we can leverage the CMAKE_SKIP_LINTING variable. The
#      ../dependencies.cmake file demonstrates how this works

# set `<out>` to 1 or 0 to indicate if `${<variable>}` refers to the value of a
# non-cache variable
#
# Context:
# - `$CACHE{<variable>}` expands to the value of the cache variable named
#   `<variable>` if it exists or to an empty string
# - `${<varname>}` to a value determined by the following order of preference:
#   1. the value of the variable named `<variable>`, if it exists
#   2. the value of the cache variable named `<variable>`, if it exists
#   3. an empty string
function(_has_noncache_definition variable out)
  if (NOT DEFINED ${variable})
    set(${out} 0 PARENT_SCOPE)
  elseif(DEFINED ${variable} AND NOT DEFINED CACHE{${variable}})
    set(${out} 1 PARENT_SCOPE)
  elseif(NOT ("${${variable}}" STREQUAL "$CACHE{${variable}}"))
    set(${out} 1 PARENT_SCOPE)
  else()
    # we know that a cache variable named <variable> definitely exists and a
    # non-cache variable named <variable> only exists if it stores exactly
    # the same value as the cache variable
    set(original_val $CACHE{${variable}})

    # we use set_property to overwrite $CACHE{${variable}} rather than set(...)
    # because the latter requires us to specify the TYPE and HELPSTRING
    if (NOT ("${original_val}" STREQUAL ""))
      set_property(CACHE ${variable} PROPERTY VALUE "")
    else()  # use a different dummy value
      set_property(CACHE ${variable} PROPERTY VALUE "<DUMMY>")
    endif()

    # now we do the comparison
    if ("${${variable}}" STREQUAL "${original_val}")
      set(${out} 1 PARENT_SCOPE)
    else()
      set(${out} 0 PARENT_SCOPE)
    endif()

    # restore the value of the cache value...
    set_property(CACHE ${variable} PROPERTY VALUE "${original_val}")
  endif()
endfunction()

# for each standard linter variable name (e.g. CMAKE_CXX_CLANG_TIDY or
# CMAKE_C_INCLUDE_WHAT_YOU_USE) that doesn't expand to an empty string,
# - create a non-CACHE variable with the same name that holds an empty string
#   (if a non-CACHE variable with that name exists, it is overridden)
# - AND save the information needed to undo this operation in a collection
#   of variables named `<prefix>_<n>` where `<n>` is an integer from 0 through
#   the value stored in `<prefix>_<COUNT>`
#
# This function is undone by passing the same `prefix` arg to the
# `LinterVar_Restore` command
#
# Implementation Notes:
# - we create separate variables for each overridden variable (rather than a
#   list) because the "standard linter variables" already hold a list (i.e.
#   semi-colons delimit extra arguments passed to the linter)
# - we choose to overwrite variables, rather than simply unset them because
#   restoring cache variables gets messy. Specifically, we would need to:
#   - track and restore the ADVANCED, HELPSTRING, STRINGS, and TYPE properties
#     of each cache variable and make sure we properly restore them.
#   - technically, we can't perfectly do this for a TYPE value of INTERNAL
#     (but I think that's ok)
#   - plus, we would have to ignore the cache variable's MODIFIED property that
#     CMake internally manages for tracking interactive modification (it's not
#     totally clear if there would be any ramifications for doing this would be)
function(LinterVar_TmpOverride prefix)
  set(counter 0)
  set(linter_list "CPPLINT;CLANG_TIDY;CPPCHECK;ICSTAT;INCLUDE_WHAT_YOU_USE")
  foreach(LINTER IN LISTS linter_list)
    foreach(LANG IN ITEMS "C" "CXX")
      set(varname "CMAKE_${LANG}_${LINTER}")
      if ((NOT DEFINED "${varname}") OR ("${${varname}}" STREQUAL ""))
        continue()
      endif()

      _has_noncache_definition(${varname} has_noncache_def)
      if ("${has_noncache_def}")
        set(old_nameval_pair "${varname}:${${varname}}")
      else()
        set(old_nameval_pair "${varname}:")
      endif()

      set(${prefix}_${counter} "${old_nameval_pair}" PARENT_SCOPE)
      set(${varname} "" PARENT_SCOPE)

      math(EXPR counter "${counter} + 1")
    endforeach()
  endforeach()

  set("${prefix}_COUNT" ${counter} PARENT_SCOPE)
endfunction()

# reverts all modifications made an earlier call to LinterVar_TmpOverride (the
# `prefix` arguments must match)
function(LinterVar_Restore prefix)
  if (NOT DEFINED "${prefix}_COUNT")
    message(FATAL_ERROR "no variable named: ${prefix}_COUNT")
  endif()

  math(EXPR loop_end "${${prefix}_COUNT} - 1")
  unset(${prefix}_COUNT PARENT_SCOPE)
  foreach(counter RANGE ${loop_end})  # upper bound of CMake range is inclusive
    # nameval_pair holds a string of the form "<varname>:<varval>""
    set(nameval_pair "${${prefix}_${counter}}")
    unset(${prefix}_${counter} PARENT_SCOPE)

    string(FIND "${nameval_pair}" ":" sep_pos)
    string(LENGTH "${nameval_pair}" len)

    string(SUBSTRING "${nameval_pair}" 0 ${sep_pos} variable_name)
    math(EXPR sep_pos "${sep_pos} + 1")

    if ("${sep_pos}" EQUAL "${len}")
      unset(${variable_name} PARENT_SCOPE)
    else()
      string(SUBSTRING "${nameval_pair}" ${sep_pos} ${len} variable_val)
      set(${variable_name} ${variable_val} PARENT_SCOPE)
    endif()
  endforeach()
endfunction()


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

