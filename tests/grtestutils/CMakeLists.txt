# Purpose: define "libraries" to aid with the testing of Grackle
# ==============================================================

# first, define an interface "library" (i.e. nothing is compiled) for
# book-keeping purposes. dependents get "linked" against this library in order
# to be able to access grackle's internal functionality

add_library(_grackle_internals INTERFACE)
target_link_libraries(_grackle_internals INTERFACE Grackle::Grackle)

# short-term macro definitions (both will become unnecessary in the near future)
target_compile_definitions(_grackle_internals
  # currently required for functions declared by fortran_func_wrappers.hpp
  INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>"
  # suppresses warnings about C internal headers using our deprecated public
  # headers (the files should include grackle.h instead). We're currently
  # holding off on properly fixing this to minimize conflicts with pending PRs
  # on the newchem-cpp branch
  INTERFACE GRIMPL_PUBLIC_INCLUDE=1
)

# this makes it possible to write an include-directive of a header from
# src/clib, by writing something like
#   #include "{name}.hpp"
# where {name} is replaced by the name of the header
# -> frankly, I don't love this because it's not immediately clear that we are
#    including a private header. The alternative is to require the paths
#    in the include directives to include the name of the directory holding
#    the implementation files (e.g. "clib/{header}.hpp")
target_include_directories(_grackle_internals
  INTERFACE ${PROJECT_SOURCE_DIR}/src/clib
)

# grtestutils: testing utility library
# ====================================
# - this is an internal library that defines reusable testing utilities
# - Frankly, I don't love the name. Originally I called it grtest, but I'm a
#   a little concerned that may look a little too much like gtest (i.e. the
#   shorthand that googletest uses).

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
#   functionality can be introduced. The organization strategy ensures that
#   it will be straightforward to provide this hypothetical functionality to
#   the python bindings (and any C++ testing code) without making the python
#   bindings link against googletest.

add_library(grtestutils
  # files outside of the googletest subdirectory
  # -> these shouldn't include headers from the googletest subdirectory
  cmd.hpp cmd.cpp
  iterator_adaptor.hpp
  os.hpp os.cpp
  utils.hpp utils.cpp

  harness/grackle_ctx_pack.cpp harness/grackle_ctx_pack.hpp
  harness/param.cpp harness/param.hpp
  harness/preset.cpp harness/preset.hpp
  harness/status.cpp harness/status.hpp

  # files in the googletest subdirectory (reminder: should just be headers!)
  googletest/check_allclose.hpp
  googletest/fixtures.hpp
)

target_link_libraries(grtestutils
  PRIVATE _grackle_internals
  PUBLIC Grackle::Grackle
  # to express the transitive dependency on googletest (it's transitive since
  # its only referenced by header files), we should theoretically write
  #   INTERFACE GTest::GTest
  # we choose not to do this since all test code pulls in this dependency by
  # explicitly listing `testdeps` as a dependency (defined below)
)

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


# 2. testdeps: bundles dependencies used by tests
# ===============================================
# - this is an interface library (i.e. nothing is compiled) that ONLY exists
#   to aggregate dependencies, include directories and compile definitions
# - all tests should depend on this target

add_library(testdeps INTERFACE)
target_link_libraries(testdeps
  INTERFACE Grackle::Grackle grtestutils _grackle_internals
            GTest::gtest_main GTest::gmock
)
# indicates that all consumers are written in C++17 or newer
target_compile_features(testdeps INTERFACE cxx_std_17)


