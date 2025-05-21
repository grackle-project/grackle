# Do 2 things before we declare the project
# 1. record the languages enabled by any parent projects
get_property(PRE_ENABLED_LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)

# 2. load and and set the version information
set(grackle_vfile ${CMAKE_CURRENT_SOURCE_DIR}/VERSION)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${grackle_vfile})

file(READ ${grackle_vfile} grackle_vfile_content)
if (grackle_vfile_content MATCHES
    "^([0-9]+(\\.[0-9_]+)*(-[-+A-Za-z0-9_]+)?)[ \t\r\n]*$")
  set(_GRACKLE_FULL_VERSION "${CMAKE_MATCH_1}")
  string(REGEX REPLACE "(-[-+A-Za-z0-9_]*[-+A-Za-z_][-+A-Za-z0-9_]*)" ""
    _GRACKLE_VERSION_NO_DEVTAG "${_GRACKLE_FULL_VERSION}")
else()
  message(FATAL_ERROR "Invalid version number in ${grackle_vfile}")
endif()


# Preamble
# --------
cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0077 NEW)
cmake_policy(SET CMP0082 NEW)

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
  cmake_policy(SET CMP0135 NEW) # makes automatic dependency downloads more
                                # robust (& addresses noisy warnings)
endif()

project(Grackle VERSION ${_GRACKLE_VERSION_NO_DEVTAG} LANGUAGES C CXX Fortran)

# Specify that `./cmake` contains cmake modules defining useful functionality
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# determine whether we are building a stand-alone program, or if grackle has
# been embedded within another project
set(GRACKLE_IS_TOP_LEVEL ON)
if (NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(GRACKLE_IS_TOP_LEVEL OFF)
endif()

# Project-Wide Configuration
# --------------------------

# perform some basic checks
if (NOT GRACKLE_IS_TOP_LEVEL)
  if ( NOT ("Fortran" IN_LIST PRE_ENABLED_LANGUAGES) )
    # technically, this is only an issue when building Grackle as a static
    # library, but we raise an error anyways for a consistent experience
    #
    # based on the answers to the following post, we will reassess the approach
    # we take in this scenario
    #   https://discourse.cmake.org/t/conventions-for-linking-implicit-dependencies-of-an-embedded-multi-language-static-library/11073
    message(FATAL_ERROR
      "The top-level project in which the Grackle-build is embedded doesn't"
      "list Fortran as an enabled language. This will cause linking issues"
      "when Grackle is built as a static library (like the current build)."
      ""
      "To fix this, either add Fortran to the language list specified in"
      "the top-level project command OR call enable_language(Fortran)"
      "before embedding Grackle (see the latest CMake docs for restrictions"
      "about where enable_language may be called")
  endif()

elseif (GRACKLE_IS_TOP_LEVEL)
  # the following 2 checks are only enforced when Grackle is built standalone
  # (we let people do whatever when it's embedded for maximum flexibility)

  if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(FATAL_ERROR "Please don't perform an in-source build.")
  endif()

  if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "no build type detected, defaulting to Release build-type")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  endif()
endif()

option(GRACKLE_USE_DOUBLE "compile the code with double precision" ON)
option(GRACKLE_DEFINE_INSTALL_RULES
  "enable Grackle's installation rules" ${GRACKLE_IS_TOP_LEVEL})
option(GRACKLE_BUILD_TESTS "build Grackle tests driven by CTest" OFF)

if (GRACKLE_BUILD_TESTS)
  set(_GRACKLE_EXAMPLES_DEFAULT ON)
else()
  set(_GRACKLE_EXAMPLES_DEFAULT ${GRACKLE_IS_TOP_LEVEL})
endif()
option(GRACKLE_EXAMPLES
  "compile grackle examples" ${_GRACKLE_EXAMPLES_DEFAULT})

# configure custom machine-specific optimization-flags (specified in host-config
# files). I'm not really sure we should support this when Grackle isn't the
# top-level project (hence the warning...)

if (NOT DEFINED GRACKLE_OPTIMIZATION_FLIST)
  set(GRACKLE_OPTIMIZATION_FLIST "${GRACKLE_OPTIMIZATION_FLIST_INIT}"
    CACHE STRING "\
machine-specific semicolon delimited list of optimization flags. These are \
only in RELEASE or RELWITHDEBINFO build-configurations. See the installation \
documentation for notes about \"shell-quoting\" groups of options. NOTE: We \
may not use these when Grackle isn't the top-level project (i.e. if it is \
embedded inside another cmake project).")
endif()

if((NOT GRACKLE_IS_TOP_LEVEL) AND
   (NOT "${GRACKLE_OPTIMIZATION_FLIST}" STREQUAL ""))
  message(WARNING
    "It's unclear whether we should support machine-specific optimization "
    "flags, when Grackle isn't the top-level project (i.e. it is embedded in "
    "the build of another project). In this scenario, it might be better to "
    "encourage the top-level project to set the ``CMAKE_<LANG>_FLAGS`` "
    "variable or use the ``add_compile_options`` command. Be aware, we may "
    "alter support for flags specified through GRACKLE_OPTIMIZATION_FLIST (& "
    "related variables) in this scenario, at any time (without any warning)")
endif()


# handle some organization for generated files (mostly for our own convenience)
include(GNUInstallDirs)
if (NOT GRACKLE_IS_TOP_LEVEL)
  set(stageDir ${CMAKE_CURRENT_BINARY_DIR})
  # the location where the cmake export files go to export from build-tree
  set(GRACKLE_BUILD_EXPORT_PREFIX_PATH ${stageDir}) # <-- may not get used

elseif(DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY OR
       DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY OR
       DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

  message(FATAL_ERROR "unexpected scenario!")

else()

  set(stageDir ${CMAKE_CURRENT_BINARY_DIR}/grackle)

  # Add an explanatory direcotory
  file(WRITE ${stageDir}/README.md [==[
This directory only exists for the convenience of the developers of Grackle.
Files may be moved around (or cease to exist). Please do not rely upon
the structure of this directory
]==])

  # these intentionally are not CACHE variables
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${stageDir}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${stageDir}/${CMAKE_INSTALL_LIBDIR})

  # don't currently need the following since grackle doesn't ship an executable
  #set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${stageDir}/${CMAKE_INSTALL_BINDIR})

  # the location where the export files go to export from build-tree
  set(GRACKLE_BUILD_EXPORT_PREFIX_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
endif()

# create a directory where we will store auto-generated generated public headers
set(GRACKLE_GENRATED_PUBLIC_HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/generated_public_headers)

# Dependencies
# ------------
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

option(GRACKLE_USE_OPENMP "Use OpenMP" OFF)
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
add_library(toolchain::m INTERFACE IMPORTED)
if (UNIX AND NOT CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set_target_properties(toolchain::m PROPERTIES IMPORTED_LIBNAME "m")
endif()

# Main build targets
# ------------------
add_subdirectory(src/clib)

# declare build-recipies for examples
if (GRACKLE_EXAMPLES)
  add_subdirectory(src/example)
endif()

# Tests
# -----

# if desired, setup Grackle's tests that use CTest/GoogleTest (this has no
# impact on the tests run via pytest, those are handled outside of cmake)
if(GRACKLE_BUILD_TESTS)
  # setup GoogleTest framework (search for an installed copy on the machine
  # OR download it & add its compilation to the rest of the build procedure)
  include(setup_GTest)
  # load CMake's built-in module that defines functions to registering unit
  # tests implemented with GoogleTest into CMake's test-runner, CTest
  include(GoogleTest)
  enable_testing() # enable the test-runner, CTest
  add_subdirectory(tests) # declare the recipes for each test
endif()

# Packaging
# ---------
# (probably not needed if not root project)
if(GRACKLE_DEFINE_INSTALL_RULES AND NOT CMAKE_SKIP_INSTALL_RULES)
  include(cmake/installation_rules.cmake)
endif()
