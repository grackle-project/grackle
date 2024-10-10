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

set(_GRDATA_CACHE_PREFIX "_GRACKLEGRDATAPRIVATE_")

# create the grdata program
#
# the program's resulting location is controlled by the global
# CMAKE_RUNTIME_OUTPUT_DIRECTORY variable, if it is set. Otherwise, we put it
# in the CURRENT_BINARY_DIR.
#
# Arguments
# ---------
#
# GRACKLE_VERSION: the full grackle version number of the grackle version that
#     the program is associated with (not just the truncated number understood
#     by cmake)
# PATH_OUTVAR: specifies the name of the variable where the path to the
#     resulting program is stored
function(create_grdata_program)

  set(options)
  set(oneValueArgs DESTINATION GRACKLE_VERSION PATH_OUTVAR)
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

  set("${CREATE_GRDATA_PATH_OUTVAR}" "${output_path}" PARENT_SCOPE)

  set(_GRACKLEGRDATAPRIVATE_TOOL_PATH "${output_path}" CACHE INTERNAL
      "cached location of the grdata tool")
endfunction()


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

# helper function to create an export-file (that will be installed) and is used
# to provide some basic properties about the grdata tool
#  EXPORT_FILE_DESTINATION_DIR specifies directory where the export-file will
#    be installed relative to the root-install path
#  TOOL_RELATIVE_INSTALL_PATH where the export-file will be
#    installed relative to the root-install path
function(_grdata_get_export_file_contents EXPORT_FILE_DESTINATION_DIR TOOL_RELATIVE_INSTALL_PATH outVar)
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

  # now, we will sanitize EXPORT_FILE_DESTINATION_DIR (and remove any trailing
  # slashes)

  string(CONFIGURE [=[
# Autogenerated file that stores the location of the grdata tool
# -> this is directly analogous to a file that would be defined with cmake's
#    install(EXPORT ...) command.
# -> since the grdata tool is a weird sort of pseudo target (i.e. it isn't
#    compiled), we store the path to the file
# -> like `install(EXPORT ...)` (and in contrast to `export(EXPORT ...)`),
#    we use a relative path to the grdata tool
set(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}/@REL_PATH_TO_PREFIX@")
set(_GRACKLE_GRDATA_TOOL_PATH "${_IMPORT_PREFIX}/@TOOL_RELATIVE_INSTALL_PATH@")
unset(_IMPORT_PREFIX)
]=] contents @ONLY
)
  set("${outVar}" "${contents}" PARENT_SCOPE)

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

  # get the kwargs for the current mode
  if("INSTALL_TOOL" STREQUAL "${mode}")
    set(_GRDATA_IEH_oneValueArgs DESTINATION COMPONENT)
  elseif("INSTALL_EXPORT" STREQUAL "${mode}")
    set(_GRDATA_IEH_oneValueArgs DESTINATION FILE TMPDIR)
  elseif("EXPORT_EXPORT" STREQUAL "${mode}")
    set(_GRDATA_IEH_oneValueArgs FILE)
  else()
    message(FATAL_ERROR
      "${_GRDATA_IEH_name} command invoked with unexpected mode: \"${mode}\""
    )
  endif()

  # parse the arguments
  cmake_parse_arguments(_GRDATA_IEH "" "${_GRDATA_IEH_oneValueArgs}" ""
                        ${ARGN})

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

  if (NOT DEFINED _GRACKLEGRDATAPRIVATE_TOOL_PATH)
    message(FATAL_ERROR
      "${_GRDATA_IEH_name}(${mode}) can only be called AFTER a call to the "
      "create_grdata_program command.")
  endif()

  # now, actually complete the command
  if ("INSTALL_TOOL" STREQUAL "${mode}")
    install(PROGRAMS ${_GRACKLEGRDATAPRIVATE_TOOL_PATH}
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

    _grdata_get_export_file_contents(
      ${_GRDATA_IEH_DESTINATION}
      ${_GRACKLEGRDATAPRIVATE_RELATIVE_INSTALL_PATH}
      "_GRDATA_IEH_EXPORT_CONTENTS"
    )

    file(WRITE "${_GRDATA_IEH_TMPDIR}/${_GRDATA_IEH_FILE}"
      "${_GRDATA_IEH_EXPORT_CONTENTS}")
    install(FILES
      "${_GRDATA_IEH_TMPDIR}/${_GRDATA_IEH_FILE}"
      DESTINATION "${_GRDATA_IEH_DESTINATION}"
    )

  elseif ("EXPORT_EXPORT" STREQUAL "${mode}")
    string(CONFIGURE [=[
# Autogenerated file that stores the location of the grdata tool
# -> this is directly analogous to a file that would be defined with cmake's
#    export(EXPORT ...) command.
# -> since the grdata tool is a weird sort of pseudo target (i.e. it isn't
#    compiled), we store the path to the file
# -> like `export(EXPORT ...)` (and in contrast to `install(EXPORT ...)`), 
#    we use an absolute path

set(_GRACKLE_GRDATA_TOOL_PATH "@_GRACKLEGRDATAPRIVATE_TOOL_PATH@")
"]=] _GRDATA_IEH_EXPORT_CONTENTS @ONLY)

  file(WRITE "${_GRDATA_IEH_FILE}" "${_GRDATA_IEH_EXPORT_CONTENTS}")

  else()
    message(FATAL_ERROR
      "something went horribly wrong within ${_GRDATA_IEH_name}(${mode})")
  endif()

endmacro()
