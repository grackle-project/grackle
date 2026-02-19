
if (GRACKLE_IS_TOP_LEVEL AND
    (CMAKE_MINIMUM_REQUIRED_VERSION VERSION_GREATER_EQUAL "3.24"))
  message(AUTHOR_WARNING
    "Reminder: now that Grackle's minimum required CMake version >=3.24:\n"
    "- replace each call to a command named `GrBackport_FetchContent_<NAME>`\n"
       "with a call to the corresponding command called `FetchContent_<NAME>`\n"
    "- replace `include(Backport_FetchContent)` with `include(FetchContent)`\n"
    "- delete `Backport_FetchContent.cmake`"
  )
elseif(GRACKLE_IS_TOP_LEVEL AND
    (CMAKE_MINIMUM_REQUIRED_VERSION VERSION_GREATER_EQUAL "4.2"))
  message(FATAL_ERROR
    "Reminder: now that Grackle's minimum required CMake version >=4.2, "
    "delete `LinterVar_TmpOverride` & `LinterVar_Restore`"
  )
endif()

# load drop-in wrappers for FetchContent_Declare & FetchContent_MakeAvailable
# that backport support for FIND_PACKAGE_ARGS kwarg
include(Backport_FetchContent)
# load LinterVar_<action> commands
include(LinterHandling)

# modify standard variables to instruct CMake not to run static-analyis linters
# on any embedded build-targets created in calls to FetchContent_MakeAvailable
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "4.2")
  _has_noncache_definition(CMAKE_SKIP_LINTING _GRACKLE_PREDEFINED_SKIP_LINTING)
  if (_GRACKLE_PREDEFINED_SKIP_LINTING)
    set(_GRACKLE_ORIG_NONCACHE_SKIP_LINTING_VAR ${CMAKE_SKIP_LINTING})
  endif()
  set(CMAKE_SKIP_LINTING OFF)
else()
  LinterVar_TmpOverride(_GRACKLE_OLD_LINTER_VARS)
endif()

# Handle external dependencies that will must be downloaded and built from
# source if pre-built copies can't be found

# NOTE: it is idiomatic to use FetchContent with a git hash rather than the
# name of a tag or branch (since the latter 2 can freely change). If we use
# the name of a tag/branch, then CMake will query GitHub on every
# subsequent build to check whether the tag/branch is still the same

if (GRACKLE_BUILD_TESTS)  # deal with testing dependencies

  # the only testing dependency is googletest.
  message(STATUS
    "If googletest can't be found, it will be downloaded & compiled"
  )

  GrBackport_FetchContent_Declare(
    googletest # the url contains the hash for v1.15.2
    URL https://github.com/google/googletest/archive/b514bdc898e2951020cbdca1304b75f5950d1f59.zip
    FIND_PACKAGE_ARGS NAMES GTest
  )

  # In case GoogleTest's is downloaded & compiled as part of the build:
  # -> Tell its build-system not to define installation rules (since we only use
  #    it to run tests from the build-directory)
  set(INSTALL_GTEST OFF)
  # -> For Windows: Prevent overriding parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

  # actually fetch & configure googletest (if it wasn't found)
  GrBackport_FetchContent_MakeAvailable(googletest)

endif() # GRACKLE_BUILD_TESTS

# restore CMake linting variables controlling linting
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL "4.2")
  if (DEFINED _GRACKLE_ORIG_NONCACHE_SKIP_LINTING_VAR)
    set(CMAKE_SKIP_LINTING ${_GRACKLE_ORIG_NONCACHE_SKIP_LINTING_VAR})
  endif()
else()
  LinterVar_Restore(_GRACKLE_OLD_LINTER_VARS)
endif()

# find all of the other dependencies
# -> we expect the caller of the build to make these available to us in one
#    form or another

find_package(HDF5 REQUIRED COMPONENTS C)
add_library(GRACKLE_HDF5_C INTERFACE IMPORTED)
target_link_libraries(GRACKLE_HDF5_C INTERFACE ${HDF5_C_LIBRARIES})
target_include_directories(GRACKLE_HDF5_C INTERFACE ${HDF5_C_INCLUDE_DIRS})
target_compile_definitions(GRACKLE_HDF5_C INTERFACE ${HDF5_C_DEFINITIONS})
if (HDF5_VERSION VERSION_LESS "1.6")
  message(FATAL_ERROR "HDF5 version 1.6 or newer is required")
elseif(HDF5_VERSION VERSION_GREATER "1.6")
  target_compile_definitions(GRACKLE_HDF5_C INTERFACE -DH5_USE_16_API)
endif()

if (GRACKLE_USE_OPENMP)
  if (CMAKE_GENERATOR STREQUAL "Ninja")
    message(WARNING
      "using Ninja with GRACKLE_USE_OPENMP=ON may cause compilation problems."
      "The issues manifest as error with finding the \"omp_lib.h\" header "
      "that is conditionally included by some Fortran source files when using "
      "CMake. The quick fix is use Makefiles. See the docs for more info"
    )
  endif()

  if(GRACKLE_EXAMPLES)
    set(_GRACKLE_OMP_COMPONENTS C Fortran CXX)
  else()
    set(_GRACKLE_OMP_COMPONENTS C Fortran)
  endif()
  find_package(OpenMP REQUIRED COMPONENTS ${_GRACKLE_OMP_COMPONENTS})
endif()

# define target to link the math functions of the C standard library
# (i.e. the -lm flag). This is commonly needed on unix-like platforms
# -> For platforms that don't need libm, this target acts as a dummy
#    placeholder (that does nothing)
# -> The -lm flag should NOT be used on MacOS (while CMake is smart enough to
#    not pass it to the linker, it will mess with exporting linker flags)
#
# NOTE: when we start using C++ in the core grackle library, we can remove
# everything related to the toolchain::m variable (since the C++ runtime
# library is ALWAYS linked to the math functions)
add_library(toolchain::m INTERFACE IMPORTED)
if (UNIX AND NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set_target_properties(toolchain::m PROPERTIES IMPORTED_LIBNAME "m")
endif()
