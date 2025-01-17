# This is a cmake module that defines the logic for creating the grdata
# "program"
#
# In a vacuum, it would be more idiomatic to take the logic in this file and
# directly embed it in a CMakeLists.txt file located in close proximity to the
# template file.
#
# However the situation is complicated by the following 2 factors:
# 1. grdata.py is used both as a template file and as part of the pygrackle
#    package. The grdata.py file is written in such away that it can be used
#    without as part of pygrackle without any sort of variable substitution
#
# 2. Moreover, the pygrackle package is built with the scikit-build-core
#    backend, which requires a set of CMake files. Since there is no reason to
#    ever directly execute this file's logic as part of building the python
#    package it makes even more sense to try to keep the logic separate


# load the file_registry information for the current version of grackle into a
# string held by the variable specified by the outvar argument
function(load_file_registry_string UPDATE_CONFIGURE_DEPENDS outvar)
  set(path
    "${PROJECT_SOURCE_DIR}/src/python/pygrackle/file_registry/file_registry.txt"
  )
  if ("${UPDATE_CONFIGURE_DEPENDS}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}")
  endif()
  file(READ ${path} contents)
  set("${outvar}" "${contents}" PARENT_SCOPE)
endfunction()

# create the grdata program (and represent it as an executable target)
#
# To install the program, and properly expose export information, see the
# grdata_install_export_helper_ command
#
# the program's resulting location is controlled by the global
# CMAKE_RUNTIME_OUTPUT_DIRECTORY variable, if it is set. Otherwise, we put it
# in the CURRENT_BINARY_DIR.
#
# Arguments
# ---------
# GRACKLE_VERSION: the full grackle version number of the grackle version that
#                  the program is associated with (not just the truncated
#                  number understood by cmake)
# TARGET_NAME: specifies the name of the target that is created to represent
#              the program.
#
# Notes
# -----
# The grdata program can be very useful for downstream codes for testing
# purposes (and for the hypothetical scenario where we support downloading
# precompiled copies of Grackle). Due to the nature of the grdata program,
# providing the program details (i.e. namely providing the path to it) to
# downstream projects is not straight-forward.
# - In more detail, the program is a little weird in a CMake context since we
#   want people to essentially think of it as a generic command line program 
#   even though it isn't technically a compiled program (people shouldn't care
#   that it is actually a portable, executable python script).
# - I was originally hesitant to declare it as an IMPORTED executable. I was
#   primarily concerned about the scenario where other CMake-built projects
#   might consume Grackle by embedding it as a part of the build. Since the
#   documentation explains that IMPORTED executable machinery is intended to
#   represent machinery from outside the CMake, I was concerned that we could
#   encounter some unforseen side-effects.
# - Prior to CMake 3.0, we might have instead provided the tool's path in a
#   package variable (like GRACKLE_GRDATA_TOOL_PATH). While the details of how
#   we achieve this in modern CMake differ, we can still provide the path in a
#   manner similar to variable-access
#   - it is still possible to provide package variables, but the *Professional
#     CMake* book, by Craig Scott (a primary CMake developer), makes is clear
#     we should to avoid package-variables to make Grackle easily consumable.
#     In the book's 18th edition this advice is in section 40.4.
#   - essentially, we would need to take extra steps for every package variable
#     that we introduce to support downstream projects employ newer dependency
#     management machinery introduced in 2022 (cmake 3.24).
#   - the advice is to make this information accessible through a function
#     (that returns values for known keys) or as properties of a target. We
#     experimented with the function approach -- see commits just before this
#     documentation was written. It requires a somewhat involved solution
#     to completely avoid global/package variables.
#
# After giving this some more thought, it became clear that the IMPORTED
# executable approach is superior for 2 reasons:
# 1. Even if there is some unforseen side-effect of the IMPORTED executable, we
#    wiil be no worse-off than the variable-like approach.
#    - The worst imaginable side-effect is that somebody consuming Grackle
#      in an embedded manner, might find that some of CMake's source-file
#      dependency magic doesn't work right after updating the files that
#      compose the grdata tool.
#    - This scenario seems extremely pathological (and should probably never
#      arise), since it could only occur if people alter their source files
#      based on the output of the grdata tool.
#    - Regardless of the practicality, the variable-like approach definitely
#      wouldn't offer any benefits here (the same issues would still occur)
# 2. The use of targets is generally more idiomatic than variables. While we
#    still require some custom logic to get everything to work right, we need
#    less of it than the alternative. Moreover, the custom-code will be much
#    more directly analogous to standard cmake logic (making it easier to
#    understand)
function(create_grdata_program)

  set(options)
  set(oneValueArgs DESTINATION GRACKLE_VERSION TARGET_NAME)
  set(multiValueArgs)

  cmake_parse_arguments(PARSE_ARGV 0
    CREATE_GRDATA "${options}" "${oneValueArgs}" "${multiValueArgs}")

  # some basic error-handling
  set(_funcname "create_grdata_program")
  if (DEFINED CREATE_GRDATA_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "${_funcname} recieved invalid arguments: "
      "\"${CREATE_GRDATA_UNPARSED_ARGUMENTS}\"")
  elseif (DEFINED CREATE_GRDATA_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR
      "${_funcname} received the ${CREATE_GRDATA_KEYWORDS_MISSING_VALUES} "
      "keyword(s) without any associated arguments.")
  endif()

  if (DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(output_path "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/grdata")
  else()
    set(output_path "${CMAKE_CURRENT_BINARY_DIR}/grdata")
  endif()

  load_file_registry_string(TRUE "_GRDATA_FILE_REGISTRY_CONTENTS")
  set(_GRDATA_GRACKLE_VERSION "${CREATE_GRDATA_GRACKLE_VERSION}")
  configure_file(
    "${PROJECT_SOURCE_DIR}/src/python/pygrackle/utilities/grdata.py"
    "${output_path}"
    @ONLY
  )

  if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.19")
    file(CHMOD ${output_path} FILE_PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ
      WORLD_EXECUTE
    )
  else()
    execute_process(COMMAND chmod a+rx ${output_path})
  endif()

  set(_GRACKLEGRDATAPRIVATE_TARGETNAME "${CREATE_GRDATA_TARGET_NAME}"
    CACHE INTERNAL "name of the target representing the grdata tool")

  add_executable(${_GRACKLEGRDATAPRIVATE_TARGETNAME} IMPORTED GLOBAL)
  set_target_properties(${_GRACKLEGRDATAPRIVATE_TARGETNAME} PROPERTIES
    IMPORTED_LOCATION "${output_path}"
  )

endfunction()

# Helper Functions used to define export files
# --------------------------------------------

# note that cmake normalizes paths on all platforms to use forward-slashes
function(_num_path_segments path outCount)
  set(arg_descr "the path argument, of the _num_path_segments command,")
  if (path MATCHES "^/.*")
    message(FATAL_ERROR "${arg_descr} can't start with `/`")
  elseif (path STREQUAL ".")
    message(FATAL_ERROR "${arg_descr} can't currently be `.`")
  elseif ((path MATCHES "^\./.*$") OR (path MATCHES "^.*/\./.*$"))
    message(FATAL_ERROR "${arg_descr} can't currently hold a `./` segment")
  elseif (path STREQUAL "..")
    message(FATAL_ERROR "${arg_descr} can't currently be `..`")
  elseif ((path MATCHES "^\.\./.*$") OR (path MATCHES "^.*/\.\./.*$"))
    message(FATAL_ERROR "${arg_descr} can't currently hold a `../` segment")
  elseif (path STREQUAL "")
    set(${outCount} "0" PARENT_SCOPE)
    return()
  endif()

  set(count 1)
  set(remainder "${path}")

  while(NOT (remainder MATCHES "^[^/]+/*$"))  # exit loop if no slashes or all
                                              # slashes are trailing
    math(EXPR count "${count} + 1")

    # remove trailing slash(es)
    if (remainder MATCHES "^(.*/[^/])/+$")
      set(remainder "${CMAKE_MATCH_1}")
    endif()
    get_filename_component(remainder "${remainder}" DIRECTORY)
  endwhile()

  set(${outCount} "${count}" PARENT_SCOPE)
endfunction()

# create a relocatable export-file (to be placed in the installation directory)
# that declares a target representing the grdata tool
#
# Arguments
# ---------
# EXPORT_FILE_DESTINATION_DIR specifies directory where the export-file will
#                             be installed, relative to the root-install path
# TOOL_RELATIVE_INSTALL_PATH  specifies path where the grdata tool will be
#                             installed, relative to the root-install path
# TMP_FILE_LOCATION           Where to put the file (right after we create it)
function(_grdata_write_installdir_export_file
    EXPORT_FILE_DESTINATION_DIR TOOL_RELATIVE_INSTALL_PATH TMP_FILE_LOCATION
)
  # a sanity check!
  if ((EXPORT_FILE_DESTINATION_DIR MATCHES "^/.*") OR
      (TOOL_RELATIVE_INSTALL_PATH MATCHES "^/.*"))
    message(
      FATAL_ERROR 
      "_grdata_get_export_file_contents can't handle an argument that starts "
      "with a forward slash")
  endif()

  _num_path_segments("${EXPORT_FILE_DESTINATION_DIR}" num_segments)

  if (num_segments EQUAL 0)
    set(REL_PATH_TO_PREFIX "")
  else()
    string(REPEAT "../" "${num_segments}" REL_PATH_TO_PREFIX)
  endif()

  set(template [======================[
# Autogenerated file that stores the location of the grdata tool
# -> this is directly analogous to a file that would be defined with cmake's
#    install(EXPORT ...) command.
# -> since the grdata tool is a weird sort of pseudo target (i.e. it isn't
#    compiled), we store the path to the file
# -> like `install(EXPORT ...)` (and in contrast to `export(EXPORT ...)`),
#    we use a relative path to the grdata tool

set(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}/@REL_PATH_TO_PREFIX@")

add_executable(@_GRACKLEGRDATAPRIVATE_TARGETNAME@ IMPORTED)
set_target_properties(@_GRACKLEGRDATAPRIVATE_TARGETNAME@ PROPERTIES
  IMPORTED_LOCATION "${_IMPORT_PREFIX}/@TOOL_RELATIVE_INSTALL_PATH@"
)

unset(_IMPORT_PREFIX)
]======================]
)

  string(CONFIGURE "${template}" contents @ONLY)
  file(WRITE ${TMP_FILE_LOCATION} "${contents}")

endfunction()

# helper function to create a export-file to support find_package using the
# build-directory (analogous to the of `export(EXPORT)` command)
function(_grdata_write_builddir_export_file IMMEDIATE_EXPORT_FILE_PATH)

  set(template [======================[
# Autogenerated file that stores the location of the grdata tool
# -> this is directly analogous to a file that would be defined with cmake's
#    export(EXPORT ...) command.
# -> since the grdata tool is a weird sort of pseudo target (i.e. it isn't
#    compiled), we store the path to the file
# -> like `export(EXPORT ...)` (and in contrast to `install(EXPORT ...)`), 
#    we use an absolute path

set(_GRACKLE_GRDATA_TOOL_PATH "@absolute_tool_path@")
]======================]
)
  get_target_property(absolute_tool_path
    ${_GRACKLEGRDATAPRIVATE_TARGETNAME} IMPORTED_LOCATION)

  string(CONFIGURE "${template}" contents @ONLY)
  file(WRITE ${IMMEDIATE_EXPORT_FILE_PATH} "${contents}")

endfunction()

# due to the "weird" nature of the grdata tool, (we treat it like its an
# executable even though it isn't compiled), this is a command to help with
# standard export/installation options
#
# This command operates in 3 modes:
#
# 1. In ``INSTALL_TOOL`` mode, the following command
#      ```
#      grdata_install_export_helper_(INSTALL_TOOL DESTINATION <dir>
#                                    COMPONENT <component>)
#      ```
#    acts as a substitute for the `install(TARGETS grdata-target ...)` command,
#    if the grdata tool were a typicial cmake target (i.e. that got compiled).
#    In practice, this does some bookkeeping and wraps the following command:
#      ```
#      install(PROGRAMS path/to/grdata DESTINATION <dir> COMPONENT <component>)
#      ```
#
# 2. In ``INSTALL_EXPORT`` mode, the following command
#      ```
#      grdata_install_export_helper_(INSTALL_EXPORT DESTINATION <dir>
#                                    FILE <name>.cmake TMPDIR <tmp-dir>)
#      ```
#    acts as an analog to the ``install(EXPORT ...)`` command, but with 1 minor
#    difference:
#      - the ``install(EXPORT ...)`` command doesn't obviously do anything at
#        configuration time (in practice it does create the export files and
#        stores them within <build>/CMakeFiles/Export)
#      - we explicitly generate the export files at configuration time & store
#        them within directory specified by TMPDIR
#    NOTE: the namespace used when we originally defined the target is reused
#
# 3. In ``EXPORT_EXPORT`` mode, the following command
#      ```
#      grdata_install_export_helper_(EXPORT_EXPORT FILE <filename>.cmake)
#      ```
#    acts as an analog to the ``export(EXPORT ...)`` command
#
macro(grdata_install_export_helper_ mode)
  # we only take 1-value arguments (after the mode argument)

  set(_GRDATA_IEH_name "grdata_install_export_helper_")

  # get the kwargs for the current mode (all kwargs expect a single arg)
  if("INSTALL_TOOL" STREQUAL "${mode}")
    set(_GRDATA_IEH_Args DESTINATION COMPONENT)
  elseif("INSTALL_EXPORT" STREQUAL "${mode}")
    set(_GRDATA_IEH_Args DESTINATION FILE TMPDIR)
  elseif("EXPORT_EXPORT" STREQUAL "${mode}")
    set(_GRDATA_IEH_Args FILE)
  else()
    message(FATAL_ERROR
      "${_GRDATA_IEH_name} command invoked with unexpected mode: \"${mode}\""
    )
  endif()

  # parse the arguments
  cmake_parse_arguments(_GRDATA_IEH "" "${_GRDATA_IEH_Args}" "" ${ARGN})

  # check argument validity
  if (DEFINED _GRDATA_IEH_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "${_GRDATA_IEH_name}(${mode}) recieved invalid arguments: "
      "\"${_GRDATA_IEH_UNPARSED_ARGUMENTS}\"")
  elseif (DEFINED _GRDATA_IEH_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR
      "${_GRDATA_IEH_name}(${mode}) received the "
      "${_GRDATA_IEH_KEYWORDS_MISSING_VALUES} keyword(s) without any "
      "associated arguments.")
  endif()
  foreach(_GRDATA_IEH_ARG IN ITEMS ${_GRDATA_IEH_Args})
    if (NOT DEFINED "_GRDATA_IEH_${_GRDATA_IEH_ARG}")
      message(FATAL_ERROR
        "${_GRDATA_IEH_name}(${mode}) requires the `${_GRDATA_IEH_ARG}` kwarg")
    endif()
  endforeach()

  # check a precondition
  if (NOT DEFINED _GRACKLEGRDATAPRIVATE_TARGETNAME)
    message(FATAL_ERROR
      "${_GRDATA_IEH_name}(${mode}) can only be called AFTER a call to the "
      "create_grdata_program command.")
  endif()

  # now, actually complete the command
  if ("INSTALL_TOOL" STREQUAL "${mode}")
    install(PROGRAMS $<TARGET_FILE:${_GRACKLEGRDATAPRIVATE_TARGETNAME}>
      COMPONENT ${_GRDATA_IEH_COMPONENT}
      DESTINATION ${_GRDATA_IEH_DESTINATION}
    )

    set(_GRACKLEGRDATAPRIVATE_RELATIVE_INSTALL_PATH
        "${_GRDATA_IEH_DESTINATION}/grdata" CACHE INTERNAL
        "install location of the grdata tool relative to base install path")

  elseif("INSTALL_EXPORT" STREQUAL "${mode}")
    if (NOT DEFINED _GRACKLEGRDATAPRIVATE_RELATIVE_INSTALL_PATH)
      message(FATAL_ERROR
        "${_GRDATA_IEH_name}(${mode}) can only be called AFTER a call to the "
        "${_GRDATA_IEH_name}(INSTALL_TOOL) command.")
    endif()

    # create the export file
    _grdata_write_installdir_export_file(
      # where we'll put export file during install (relative to install-prefix)
      ${_GRDATA_IEH_DESTINATION}
      # where we'll put grdata tool during install (relative to install-prefix)
      ${_GRACKLEGRDATAPRIVATE_RELATIVE_INSTALL_PATH}
      # where in the build-directly we'll put the export file immediately after
      # we create it (we will copy it from here when we install)
      "${_GRDATA_IEH_TMPDIR}/${_GRDATA_IEH_FILE}"
    )

    # define the rule to copy the export-file during installation
    install(FILES
      "${_GRDATA_IEH_TMPDIR}/${_GRDATA_IEH_FILE}"
      DESTINATION "${_GRDATA_IEH_DESTINATION}"
    )

  elseif ("EXPORT_EXPORT" STREQUAL "${mode}")
    _grdata_write_builddir_export_file("${_GRDATA_IEH_FILE}")

  else()
    message(FATAL_ERROR
      "something went horribly wrong within ${_GRDATA_IEH_name}(${mode})")
  endif()

  # need to do some cleanup (since we are in a macro):
  foreach(_GRDATA_IEH_ARG IN LISTS _GRDATA_IEH_Args)
    unset("_GRDATA_IEH_${_GRDATA_IEH_ARG}")
  endforeach()

endmacro()
