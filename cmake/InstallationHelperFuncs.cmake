# assists with language-specific implicit requirements
function(_try_find_library name_l path resultVar)
  unset(LIBSEARCH_RSLT CACHE)

  if ("${path}" STREQUAL "") 
    #message("searching for ${name_l} in default search paths")
    find_library(LIBSEARCH_RSLT NAMES ${name_l} NAMES_PER_DIR)
  else()
    #message("searching for ${name_l} in ${path}")
    find_library(LIBSEARCH_RSLT
      NAMES ${name_l} NAMES_PER_DIR PATHS "${path}" NO_DEFAULT_PATH)
  endif()
  #message("rslt: ${LIBSEARCH_RSLT}")

  if (LIBSEARCH_RSLT)
    set("${resultVar}" 1 PARENT_SCOPE)
  else()
    set("${resultVar}" 0 PARENT_SCOPE)
  endif()

  unset(LIBSEARCH_RSLT CACHE)
endfunction()

# this is used to get (non-C) language-specific implicit requirements
function(get_implicit_link_reqs lang out_libs out_linkdirs)
  set(libs ${CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES})
  list(REMOVE_DUPLICATES libs)

  # remove some common dependencies (that aren't really language-specific)
  set(nondeps
    c     # libc.so     -> this is a dependency of (almost) everything
    gcc_s # libgcc_s.so -> dependency commonly needed by underlying linker
    gcc   # libgcc.a    -> dependency commonly needed by underlying linker
  )
  list(REMOVE_ITEM libs ${nondeps})

  # remove some common dependencies that we know are unnecessarily linked
  # when compiling with gfortran on macOS (i.e. they are a little redundant)
  # -> this isn't strictly necessary, but it provide a slightly cleaner result
  if (("gfortran" IN_LIST libs) AND (CMAKE_SYSTEM_NAME STREQUAL "Darwin"))
    list(REMOVE_ITEM libs "emutls_w" "heapt_w")
  endif()

  set(${out_libs} ${libs} PARENT_SCOPE)

  set(linkdirs ${CMAKE_${lang}_IMPLICIT_LINK_DIRECTORIES})

  # figure out which (if any) of directories in linkdirs needs to be used for
  # linking the implicit required directories
  # -> For each implicitly required directory, we first check for it in the
  #    implicit link directories.
  # -> To provide unsurprising behavior, we try to mimic the conventions of a 
  #    generic ``ld`` implementation. As I understand it, they generally:
  #    - look through the search directories in the order specified before
  #      checking the standard system search paths
  #    - for each directory, they look for the dynamic-version and then they
  #      search for the static version (only after this will they move on to
  #      another directory)
  #    - this is definitely how ld works on macOS since Xcode4 (2010) and it
  #      sounds like how GNU ld works (from the man page)
  # -> We will report an error if a given library can't be found anywhere
  # -> There are 2 areas with room for improvement:
  #    1. we are playing a little fast and loose about how the "standard
  #       search path" is defined (we currently assume that is the same as
  #       what cmake checks -- that is not strictly true)
  #    2. in the case of a static library, we probably want to prefer
  #
  # NOTE: WE JUST CARE ABOUT LOCATING LIBRARIES AT LINKTIME! WE DO NOT CARE
  #       ABOUT WHERE THE DYNAMIC LOADER FINDS THE LIBRARY AT RUNTIME (THAT'S
  #       A SEPARATE ISSUE...). ANYTHING MORE SOPHISTICATED REQUIRES US TO DO
  #       THIS ANYWAYS
  # find needed directories!
  set(needed_dirs "")
  foreach(lib IN LISTS libs)
    set(found_in_dir "")

    foreach(d IN LISTS linkdirs)
      _try_find_library("${lib};lib${lib}.a" "${d}" found)
      if (${found})
        set(found_in_dir "${d}")
        break()
      endif()
    endforeach()

    if ("${found_in_dir}" STREQUAL "")
      # sanity check: confirm that we can find the directory in the standard
      # search path
      _try_find_library("${lib};lib${lib}.a" "" found)
      if (NOT "${found}")
        message(FATAL_ERROR
          "Something is wrong, the implicitly required library, ${lib}, could "
          "not be found anywhere!"
        )
      endif()
    else()
      #message(STATUS "${lib} found in ${found_in_dir}")
      list(APPEND needed_dirs "${found_in_dir}")
    endif()

  endforeach()

  # create retained_dirs -- it includes all directories from linkdirs
  # that are also listed in needed_dirs, while retaining the original order
  set(retained_dirs "")
  foreach(d IN LISTS linkdirs)
    if(("${d}" IN_LIST needed_dirs) AND (NOT ("${d}" IN_LIST retained_dirs)))
      list(APPEND retained_dirs "${d}")
    endif()
  endforeach()

  set(${out_linkdirs} ${retained_dirs} PARENT_SCOPE)
endfunction()


# checks whether the hdf5.pc file exists
function(check_h5pc_exists version_req outVar)
  set(exists "FALSE")
  find_package(PkgConfig QUIET)
  if (PKG_CONFIG_FOUND)

    set(args "--modversion;hdf5")
    if (NOT "${required_version}" STREQUAL "")
      set(args "--modversion;hdf5 = ${version_req}")
    endif()

    # we could be a little more careful about version numbers...
    # we explicitly avoid using the pkg-config functions provided by PkgConfig
    # -> those define a lot of variables we don't care about
    # -> we just want to check if an hdf5 pkg-config file exists in standard
    #    search paths
    execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} ${args}
      RESULT_VARIABLE _HDF5_PC_FILE_FOUND
      OUTPUT_VARIABLE _dummy_stdout_variable
      ERROR_VARIABLE _dummy_stderr_variable
    )

    if(_HDF5_PC_FILE_FOUND EQUAL 0)
      set(exists "TRUE")
    endif()
  endif()

  set("${outVar}" "${exists}" PARENT_SCOPE)
endfunction()

# the following function actually creates grackle's pkg-config file
#
#   configure_pkgconfig_file(DESTINATION <path/to/output.pc>
#                            [STATIC_ONLY]
#                            STATIC_LIBS static_libs_str
#                            [STATIC_REQUIRES static_requires_str]
#                            INFO_PROPERTIES info_props)
#
# -> by default, we configure the file so that it is well-suited for an
#    installation with a lone dynamic library or an installation with both a
#    dynamic and static lib (here we follow all standard conventions!)
# -> by passing the ``STATIC_ONLY`` argument, we create a file that should be
#    included in an installation that only includes a static-library
#
# ``STATIC_LIBS`` specifies linker flags that are required when linking Grackle
# as a static library to a downstream application. This omits flags that are
# also needed when linking Grackle as a shared library
#
# ``STATIC_REQUIRES`` Specifies the contents of the pkg-config
# ``Requires.private`` to designate dependencies with their own pkg-config files
# that must be explicitly linked when using Grackle as a static library.
#
# ``INFO_PROPERTIES`` Specifies a properly-formatted snippet of text for the
# output file that defines a series of (queryable) project-specific variables
# that encode useful metadata
function(configure_pkgconfig_file)
  set(options STATIC_ONLY)
  set(oneValueArgs DESTINATION STATIC_LIBS STATIC_REQUIRES INFO_PROPERTIES)
  set(multiValueArgs)

  # cmake_parse_arguments in 2 forms, but we can't use the other one if we want
  # to support passing an empty string to STATIC_REQUIRES
  # -> in detail the other form requires us to explicitly pass in each argument
  #    by "evaluating" the ARGN variable (an auto-defined list of arguments)
  # -> when STATIC_REQUIRES is passed an empty "", then the contents of ARGN
  #    might look like "...;STATIC_REQURES;;...". The problem is that when you
  #    "evaluate" a variable for argument-forwarding, CMake treats multiple
  #    contiguous semi-colons as a single semi-colon. In other words, it looks
  #    like STATIC_REQUIRES doesn't receive any arguments
  # -> the form that we currently use doesn't involve this list evaluation
  cmake_parse_arguments(PARSE_ARGV 0 CONFIGURE_PC "${options}"
                        "${oneValueArgs}" "${multiValueArgs}")

  # some basic error-handling
  set(_funcname "configure_pkgconfig_file")
  if (DEFINED CONFIGURE_PC_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR
      "${_funcname} recieved invalid arguments: "
      "\"${CONFIGURE_PC_UNPARSED_ARGUMENTS}\"")
  elseif (DEFINED CONFIGURE_PC_KEYWORDS_MISSING_VALUES)
    message(FATAL_ERROR
      "${_funcname} received the ${CONFIGURE_PC_KEYWORDS_MISSING_VALUES} "
      "keyword(s) without any associated arguments.")
  endif()

  foreach(kw IN LISTS oneValueArgs multiValueArgs)
    if (kw STREQUAL "STATIC_REQUIRES")
      # STATIC_REQUIRES is an optional kwarg. Additionally, when it's passed an
      # empty string, cmake_parse_arguments doesn't define the corresponding
      # variable (it appears like the kwarg wasn't specified at all)
      continue()
    elseif (NOT DEFINED CONFIGURE_PC_${kw})
      message(FATAL_ERROR "${_funcname} didn't receive the ${kw} argument")
    endif()
  endforeach()

  # now do the heavy lifting!
  set(_PC_INFO_PROPERTIES "${CONFIGURE_PC_INFO_PROPERTIES}")

  if (CONFIGURE_PC_STATIC_ONLY)
    # the leading space on the next line is important
    set(_PC_EXTRA_LIBS " ${CONFIGURE_PC_STATIC_LIBS}")

    if (CONFIGURE_PC_STATIC_REQUIRES STREQUAL "")
      set(_PC_REQUIRES_ENTRY "\
# Requires omitted (pkg-config file(s) weren't found for any dependency,
#                   Libs.private is modified accordingly)")
    else()
      set(_PC_REQUIRES_ENTRY "Requires: ${CONFIGURE_PC_STATIC_REQUIRES}")
    endif()
    set(_PC_LIBS_PRIVATE_ENTRY "\
# Libs.private omitted (unneeded since this file describes a static library in
#                       the absence of a shared library)")
    set(_PC_REQUIRES_PRIVATE_ENTRY "\
# Requires.private omitted (unneeded since this file describes a static library
#                           in the absence of a shared library)")

  else()
    # for this case, the hardcoded `Libs` value in the template is sufficient
    set(_PC_EXTRA_LIBS "")
    set(_PC_REQUIRES_ENTRY "\
# Requires omitted (unnecessary for install with a shared lib)")

    set(_PC_LIBS_PRIVATE_ENTRY "Libs.private:${CONFIGURE_PC_STATIC_LIBS}") 
    if (CONFIGURE_PC_STATIC_REQUIRES STREQUAL "")
      set(_PC_REQUIRES_PRIVATE_ENTRY "\
# Requires.private omitted (pkg-config file(s) weren't found for any dependency,
#                           Libs.private is modified accordingly)")
    else()
      set(_PC_REQUIRES_PRIVATE_ENTRY
        "Requires.private: ${CONFIGURE_PC_STATIC_REQUIRES}")
    endif()
  endif()

  configure_file(${PROJECT_SOURCE_DIR}/cmake/grackle.pc.in
    ${CONFIGURE_PC_DESTINATION} @ONLY)
endfunction()

