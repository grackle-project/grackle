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

# 1. grtestutils: testing utility library
# =======================================
# - this is an internal library that defines reusable testing utilities
# - Frankly, I don't love the name. Originally I called it grtest, but I'm a
#   a little concerned that may look a little too much like gtest (i.e. the
#   shorthand that googletest uses).

# Code Organization
# -----------------
# This code in this directory is organized in a slightly peculiar manner
#
# Code inside and outside of the ./googletest subdirectory are compiled as
# separate (but related) libraries:
# We compile 2 libraries:
# 1. `grtestutils_core`
#    - includes all code outside of ./googletest subdirectory
#    - provides functionality useful in several contexts (see RATIONALE)
# 2. `grtestutils_googletest`
#    - **ONLY** includes code inside of ./googletest subdirectory
#    - defines custom googletest assertions and fixtures
#
# Key difference:
# - `grtestutils_core` must not include ANY dependence on the googletest library
# - `grtestutils_googletest` includes the only code within this directory that
#   depends on googletest
# - as a consequence `grtestutils_googletest` depends upon `grtestutils_core`,
#   but `grtestutils_core` MUST be not depend upon `grtestutils_googletest`
#
# RATIONALE:
# There are multiple scenarios where we want to write C++ code using Grackle
# that isn't part of the core Grackle library:
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
# - writing test-scenarios using Grackle for (i) unit tests, (ii) system tests,
#   & (iii) benchmarks (PR #506 provides a WIP implementation of the latter).
#   Writing them in C++ rather than with python bindings:
#   - is extremely useful if you want to use tools like valgrind or
#     compiler sanitizers. While it's generally possible to use these tools in
#     a python extension module, it may require compiling CPython itself from
#     source. Moreover it can be extremely slow (these tools introduce to
#     basic C/C++ operations and that overhead will be apply to the Python
#     interpretter)
#   - is important if we want to test Grackle when compiled with various
#     backends. Historically, the python bindings have never been able to use
#     Grackle with the OpenMP backend. While we can and will address this,
#     there will additional challenges as we introduce GPU acceleration
#   - makes it much easier to profile test cases (to diagnose why code is slow)
#   - provides more specific benefits:
#     - unit tests aren't very tractable with python. To test an internal
#       module with python bindings, you would need to expose each relevant
#       function/struct in a header (while making it clear it isn't part of the
#       public API), declare the function/struct for cython, and write a python
#       wrapper. This is a lot of extra friction to make changes
#     - system tests:
#       - certain classes of tests are easier to write with C++ (e.g. setting
#         metal_density to a nullptr, varying ghost zones). While such tests
#         are possible in python, they would making and maintaining a bunch of
#         extra changes to the python bindings that are only used for testing
#       - investigating certain classes of bugs is easier with C++ tests because
#         it's easier to attach gdb/lldb & introduce breakpoints
#       - (less importantly) C++ systems tests are a lot faster
#     - for benchmarks: using C++ (rather than python bindings), facillitates
#       more robust measurements (avoiding the python interpretter's overhead).
#       At this time, PR #506 provides a WIP benchmark implementation
#
# This library is intended as a location where common functionality used in
# various scenarios can be introduced. The organization strategy ensures that
# it will be straightforward to provide this hypothetical functionality to
# the python bindings, the benchmarks and any C++ testing code without making
# the bindings or benchmarks link against googletest.

# 1a. grtestutils_core
# --------------------
add_library(grtestutils_core
  # reminder:
  # - only includes files outside of the googletest subdirectory
  # - these shouldn't reference any headers from the googletest subdirectory
  cmd.hpp cmd.cpp
  iterator_adaptor.hpp
  os.hpp os.cpp
  preset.hpp preset.cpp
  ratequery.hpp
  utils.hpp utils.cpp
  view.cpp view.hpp
)

target_link_libraries(grtestutils_core
  PRIVATE _grackle_internals
  PUBLIC Grackle::Grackle
)

target_compile_features(grtestutils_core PUBLIC cxx_std_17)

# when a file outside of this directory includes a header from within this
# library, the include directive looks something like:
#    #include "grtestutils/<path/to/name>.hpp"
target_include_directories(grtestutils_core
  INTERFACE ${PROJECT_SOURCE_DIR}/tests)

# these compile-definitions act as short-term hacks
target_compile_definitions(grtestutils_core
  # this hack helps us get the path to the directory holding data files (we can
  # remove it once we introduce automatic file management in PR 235, PR 237,
  # and PR 246)
  PRIVATE GR_DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/../../grackle_data_files/input/

  # this hack lets us use Operating-system specific functionality (Once PR #237
  # is merged, we should make use of the machinery introduced by that PR for
  # enabling/disabling os-specific features)
  PRIVATE "$<$<PLATFORM_ID:Linux,Darwin>:PLATFORM_GENERIC_UNIX>"
)

# now declare targets that depend upon googletest
if (TARGET GTest::gtest)
  # 1b. grtestutils_googletest
  # --------------------------
  # - right now this is an interface library (i.e. nothing is compiled) that is
  #   header-only.
  # - The CMake target exists to comunicate build-requirements for dependents
  #   to make use of the header files within the ./googletest subdirectory
  add_library(grtestutils_googletest
    # reminder: only includes files inside of the googletest subdirectory
  
    googletest/assertions.hpp
    googletest/check_allclose.cpp googletest/check_allclose.hpp
    googletest/fixtures.hpp
    googletest/matchers.hpp
  )

  target_link_libraries(grtestutils_googletest
    PUBLIC grtestutils_core GTest::gtest
  )

  # 2. testdeps: bundles dependencies used by tests
  # ===============================================
  # - this is an interface library (i.e. nothing is compiled) that ONLY exists
  #   to aggregate dependencies, include directories and compile definitions
  # - all tests should depend on this target
  
  add_library(testdeps INTERFACE)
  target_link_libraries(testdeps
    INTERFACE Grackle::Grackle grtestutils_googletest _grackle_internals
              GTest::gtest_main GTest::gmock
  )
  # indicates that all consumers are written in C++17 or newer
  target_compile_features(testdeps INTERFACE cxx_std_17)
endif()


