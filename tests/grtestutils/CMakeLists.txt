# Purpose: Defines 2 "libraries" to aid with the testing of Grackle
# =================================================================

# 1. testdeps: bundles dependencies used by all tests
# ===================================================
# - this is an interface library (i.e. nothing is compiled) that ONLY exists
#   to aggregate dependencies, include directories and compile definitions
# - in other words, a basic unit test (that doesn't require logic from
#   grtestutils can list this as its sole dependency)

add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# compiler defs are short-term hacks to help tests invoke internal functions
# from Grackle (both of these will become unnecessary in the near future)
target_compile_definitions(testdeps
  # currently needed to invoke the fortran function wrappers from within clib
  INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>"
  # suppresses warnings about C internal headers using our deprecated public
  # headers (the files should include grackle.h instead). We're currently
  # holding off on properly fixing this to minimize conflicts with pending PRs
  # on the newchem-cpp branch
  INTERFACE GRIMPL_PUBLIC_INCLUDE=1
)

# needed to let tests call internal functions from src/clib
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)
target_compile_features(testdeps INTERFACE cxx_std_17)


# 2. grtestutils: testing utility library
# =======================================
# - this is an internal library that defines reusable testing utilities
# - Frankly, I don't love the name. Originally I called it grtest, but I'm a
#   a little concerned that may look a little too much like gtest (i.e. the
#   shorthand that googletest uses.

# Code Organization
# -----------------
# This code in this directory is organized in a slightly peculiar manner
# - code in the main directory should not depend upon googletest
# - the googletest subdirectory is intended to define custom googletest
#   assertions and fixtures. This subdirectory should **ONLY** include
#   headers (so that, at most, this internal library has a transitive
#   dependence on googletest)
#
# RATIONALE:
# - over time, the idea has been floated to introducing C++ logic to
#   accelerate certain calculations provided by the python wrapper:
#   - PR #343 proposed creating a C++ accelerated harness for performing
#     certain kinds of semi-analytic calculations (like freefall or evolving
#     constant density)
#   - PR #370 proposed introducing a function to infer the internal energy at
#     a desired temperature (yes, the PR proposes it as a python function, but
#     it would be faster as a C++ function).
#   - Relatedly, it would also be useful to be able to directly compute
#     temperature or chemical equilibrium
# - in each case this functionality would also be extremely useful for the
#   googletest framework or for writing simple C++-only benchmarking programs
#   - this is extremely useful if you want to use tools like valgrind or
#     compiler sanitizers. While it's generally possible to use these tools in
#     a python extension module, it may require compiling CPython itself from
#     source. Moreover it can be extremely slow (these tools introduce to
#     basic C/C++ operations and that overhead will be apply to the Python
#     interpretter)
#   - this is important if we want to test grackle when compiled with various
#     backends. Historically, the python bindings have never been able to use
#     Grackle with the OpenMP backend. While we can and will address this,
#     there will additional challenges as we introduce GPU acceleration
# - this library is intended as a location where that kind of hypothetical
#   functionality can be introduced, and organization strategy ensures that
#   it will be straightforward to provide this hypothetical functionality to
#   the python bindings without it depending upon googletest (i.e. we may
#   treat the googletest subdirectory as a distinct interface library)

add_library(grtestutils
  cmd.hpp cmd.cpp
  utils.hpp utils.cpp
  os.hpp os.cpp

  # REMINDER: the googletest subdirectory should only contain headers
  googletest/check_allclose.hpp
)

# we are being a little lazy with our usage of testdeps right here
target_link_libraries(grtestutils PUBLIC testdeps)
target_compile_features(grtestutils PUBLIC cxx_std_17)

# configure the grtestutils target so that when a file outside of this
# directory includes a header from within this directory, the include
# directive looks something like:
#    #include "grtestutils/<name>.hpp"
# OR
#    #include "grtestutils/googletest/<name>.hpp
target_include_directories(grtestutils INTERFACE ${PROJECT_SOURCE_DIR}/tests)

# these compile-definitions act as short-term hacks
target_compile_definitions(grtestutils
  # this hack helps us get the path to the directory holding data files (we can
  # remove it once we introduce automatic file management in PR 235, PR 237,
  # and PR 246)
  PRIVATE GR_DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/../../grackle_data_files/input/

  # this hack lets us use Operating-system specific functionality (Once PR #237
  # is merged, we should make use of the machinery introduced by that PR for
  # enabling/disabling os-specific features)
  PRIVATE "$<$<PLATFORM_ID:Linux,Darwin>:PLATFORM_GENERIC_UNIX>"
)
