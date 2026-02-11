# This module defines the wrappers:
#      wrapper-name                      |   drop-in replacement for
#  --------------------------------------|-------------------------------------
#  GrBackport_FetchContent_Declare       | FetchContent_Declare
#  GrBackport_FetchContent_MakeAvailable | FetchContent_MakeAvailable
#
# By using the wrappers, you can call GrBackport_FetchContent_Declare with the
# FIND_PACKAGE_ARGS and then GrBackport_FetchContent_MakeAvailable does the
# "right thing"
#
# In more detail:
# - when using a cmake version before 3.24, the GrBackport_FetchContent_Declare
#   wrapper intercepts the FIND_PACKAGE_ARGS kwarg, slightly modifies behavior,
#   and the GrBackport_FetchContent_MakeAvailable wrapper does the "right thing"
# - compared to the CMake =>3.24 implementation, our versions:
#   - don't understand the OVERRIDE_FIND_PACKAGE kwarg
#   - are a little more "eager" (i.e. internal find_package calls that we do
#     in GrBackport_FetchContent_Declare should technically happen in
#     GrBackport_FetchContent_MakeAvailable backported versions.
#   - don't support dependency-provider machinery
#   - this is all "good enough" for our purposes
# - when using a cmake version >= 3.24, the wrappers just directly forwards all
#   arguments onto the canonical implementations of FetchContent_Declare &
#   FetchContent_MakeAvailable
#
# Motivation:
# -> this is extremely useful for projects that assume CMake >= 3.24 (CMake 3.24
#    made some changes to streamline and unify dependency management in CMake
#    and we now do things "properly")
# -> when we eventually bump Grackle's minimum required CMake version to 3.24
#    it will be trivial to delete this file

include_guard()

include(FetchContent)

# ASIDE: perfect function arg-forwarding is a little clunky in CMake
#   -> the way that variable expansion works produces surprising results if you
#      aren't careful. We adopt the technique recommended by one of the CMake
#      maintainers (Craig Scott) from his book "Professional CMake"
#   -> the author is very clear at the start of the book that readers are free
#      to reuse sample code without attribution or licensing
#   -> for the argument forward that we are doing, we could **probably** do
#      something simpler, but it's better to be safe than sorry

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  # the wrappers simply forward args to FetchContent

  function(GrBackport_FetchContent_Declare)
    cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
    set(quotedArgs "")
    foreach(arg IN LISTS FWD_UNPARSED_ARGUMENTS)
      string(APPEND quotedArgs " [===[${arg}]===]")
    endforeach()
    cmake_language(EVAL CODE "FetchContent_Declare(${quotedArgs})")
  endfunction()

  function(GrBackport_FetchContent_MakeAvailable)
    cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
    set(quotedArgs "")
    foreach(arg IN LISTS FWD_UNPARSED_ARGUMENTS)
      string(APPEND quotedArgs " [===[${arg}]===]")
    endforeach()
    cmake_language(EVAL CODE "FetchContent_MakeAvailable(${quotedArgs})")
  endfunction()

else()

  # we make use of _GRACKLE_Backport_FetchContent_BYPASS_DOWNLOAD to communicate
  # names to Backport_FetchContent_MakeAvailable that we have already found via
  # a call to find_package from within Backport_FetchContent_Declare
  #
  # For the uninitiated: a comparison b/t GLOBAL property & a CACHE variable
  # - similarity: both make it sense to store information at global scope
  # - difference: a CACHE variable is saved between CMake runs
  # - difference: the way you get and access data is slightly different
  #
  # I'm choosing to use a GLOBAL property in this case b/c it's plausible that
  # caching this information could cause problems (it probably would be fine, but
  # I don't want to take that chance)
  
  function(GrBackport_FetchContent_Declare content_name)

    # build up quotedNormalArgs variables (it holds args that we forward)
    # -> if we encounter FIND_PACKAGE_ARGS, create and build-up a variable
    #    called quotedFindPackageArgs with all following args (if any)
    # -> otherwise, quotedFindPackageArgs is undefined
    cmake_parse_arguments(PARSE_ARGV 1 FWD "" "" "")
    set(quotedNormalArgs "[===[${content_name}]===]")
    foreach(arg IN LISTS FWD_UNPARSED_ARGUMENTS)
      if (DEFINED quotedFindPackageArgs)
        string(APPEND quotedFindPackageArgs " [===[${arg}]===]")
      elseif("${arg}" STREQUAL "OVERRIDE_FIND_PACKAGE")
        message(FATAL_ERROR "Our backport can't handle the `${arg}` keyword")
      elseif("${arg}" STREQUAL "FIND_PACKAGE_ARGS")
        set(quotedFindPackageArgs "")
      else()
        string(APPEND quotedNormalArgs " [===[${arg}]===]")
      endif()
    endforeach()

    # store effective val of FETCHCONTENT_TRY_FIND_PACKAGE_MODE in tfpmode
    string(TOUPPER "${content_name}" UPPER_CONTENT_NAME)
    if (NOT("${FETCHCONTENT_SOURCE_DIR_${UPPER_CONTENT_NAME}}" STREQUAL ""))
      set(tfpmode NEVER)
    elseif(NOT DEFINED FETCHCONTENT_TRY_FIND_PACKAGE_MODE)
      set(tfpmode OPT_IN)
    elseif(FETCHCONTENT_TRY_FIND_PACKAGE_MODE MATCHES "^OPT_IN|ALWAYS|NEVER$")
      set(tfpmode ${FETCHCONTENT_TRY_FIND_PACKAGE_MODE})
    else()
      message(FATAL_ERRROR
          "FETCHCONTENT_TRY_FIND_PACKAGE_MODE was set to a value other than "
          "OPT_IN, ALWAYS, or NEVER")
    endif()

    # possibly modify quotedFindPackageArgs
    if ("${tfpmode}" STREQUAL "ALWAYS" AND NOT DEFINED quotedFindPackageArgs)
      set(quotedFindPackageArgs "")
    elseif("${tfpmode}" STREQUAL "NEVER" AND DEFINED quotedFindPackageArgs)
      unset(quotedFindPackageArgs)
    endif()


    if (DEFINED quotedFindPackageArgs)
      cmake_language(EVAL CODE
        "find_package(${content_name} QUIET ${quotedFindPackageArgs})")
      set(found_var_name "${content_name}_FOUND")

      if (NOT DEFINED ${found_var_name})
        message(FATAL_ERROR "sanity-check: ${found_var_name} is not defined")
      elseif("${${found_var_name}}")
        # record that GrBackport_FetchContent_MakeAvailable shouldn't download
        # data for ${content_name}
        set_property(GLOBAL APPEND
          PROPERTY _GRACKLE_Backport_FetchContent_BYPASS_DOWNLOAD
          "${content_name}"
        )
        return()
      endif()
    endif()

    cmake_language(EVAL CODE
      "FetchContent_Declare(${content_name} ${quotedNormalArgs})")
  endfunction()
  
  function(GrBackport_FetchContent_MakeAvailable)
    if ("${ARGC}" EQUAL 0)
      message(FATAL_ERROR "GrBackport_FetchContent_MakeAvailable passed no arg")
    endif()

    # store the list of each content_name in the skip_list variable that was
    # found with find_package
    get_property(skip_list GLOBAL
      PROPERTY _GRACKLE_Backport_FetchContent_BYPASS_DOWNLOAD
      )

    # build up quotedArgs with every content_name not in skip_list
    cmake_parse_arguments(PARSE_ARGV 0 FWD "" "" "")
    set(quotedArgs "")
    foreach(arg IN LISTS FWD_UNPARSED_ARGUMENTS)
      if(NOT("${arg}" IN_LIST skip_list))
        string(APPEND quotedArgs " [===[${arg}]===]")
      endif()
    endforeach()

    if (NOT ("${quotedArgs}" STREQUAL ""))
      cmake_language(EVAL CODE "FetchContent_MakeAvailable(${quotedArgs})")
    endif()
  endfunction()
  
  endif()