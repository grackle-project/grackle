# this files includes some basic machinery for enabling code-coverage
#
# there is no standard way to do this. Examples include:
# - CMake Scripts (https://git.stabletec.com/other/cmake-scripts) -- this is
#   used by hdf5
# - parthenon-hpc (https://github.com/parthenon-hpc-lab/parthenon/blob/develop/cmake/CodeCov.cmake)
# - Professional CMake by Craig Scott (a CMake maintainer) suggests creating a
#   custom build-type (alongside the standard choices Debug, Release, ...)

function(target_add_coverage_flags targetName)
  # it would probably be simpler to create a custom build-type, but I don't
  # know how easy it will be disentangle CMake stuff

  if ((NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_CXX_COMPILER_ID) OR
      (NOT CMAKE_C_COMPILER_ID STREQUAL CMAKE_Fortran_COMPILER_ID))
    # it *might* be okay, but we'll need to confirm that
    message(
      WARNING "code-coverage instrumentation may not work properly if you're "
      "using compilers from different vendors. CC: ${CMAKE_C_COMPILER_ID}, "
      "CXX: ${CMAKE_CXX_COMPILER_ID}, FC: ${CMAKE_Fortran_COMPILER_ID"
    )
  endif()

  get_property(isMultiConfig GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
  if (${isMultiConfig} OR NOT CMAKE_BUILD_TYPE MATCHES "[Dd][Ee][Bb][Uu][Gg]")
    message(
      WARNING "code-coverage information may not be meaningful if not using a "
              "\"Debug\" build-type"
    )
  endif()

  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${targetName} PRIVATE
      "--coverage" "-fprofile-abs-path"
    )
  elseif (CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang|IntelLLVM")
    target_compile_options(${targetName} PRIVATE "--coverage" )
  else()
    message(
      FATAL_ERROR "Don't know coverage flags for ${CMAKE_CXX_COMPILER_ID}"
    )
  endif()

  get_target_property(targetType ${targetName} TYPE)
  if (targetType MATCHES "STATIC_LIBRARY|OBJECT_LIBRARY|INTERFACE_LIBRARY")
    set(linker_scope INTERFACE)
  else()
    set(linker_scope PRIVATE)
  endif()
  target_link_options(${targetName} ${linker_scope} "--coverage")

endfunction()
