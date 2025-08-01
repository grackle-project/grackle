# This file is only to be used when building the pygrackle python module
# -> the python build-backend, scikit-build-core, is responsible for reading in
#    this file (it defines predefines certain variables and certain commands)
# -> scikit-build-core is configured in ${GRACKLE_REPO_ROOT}/pyproject.toml
set(GRACKLE_REPO_ROOT  ${CMAKE_CURRENT_SOURCE_DIR}/../../)

# Prologue
# --------

# we require 3.21 at the recommendation of the cython-cmake python-package
#   https://github.com/scikit-build/cython-cmake
# it's ok if this exceeds the version required for the core-grackle lib (which
# is linked to the lowest CMake version provided on common Linux distributions)
# since the scikit-build-core backend has been configured to automatically
# download a newer version of CMake, when necessary
cmake_minimum_required(VERSION 3.21)

project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

# Dependencies
# ------------

# find required python components
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
# read in cython utilities (provided by cython-cmake python-package)
include(UseCython)

# deal with dependency on Grackle and define _PYGRACKLE_CONSUME_MODE
if ("$ENV{PYGRACKLE_LEGACY_LINK}" STREQUAL "classic")
  # match old setuptools-approach for the builds with the classic build-system
  # -> explicitly match legacy behavior: link against library in local build &
  #    HOPE that the library will be installed when we need it
  # -> more idiomatic alternative: provide cmake with full library path
  set(_PYGRACKLE_CONSUME_MODE "ExternalClassic")

  add_library(Grackle::Grackle INTERFACE IMPORTED)
  target_include_directories(Grackle::Grackle
    INTERFACE ${GRACKLE_REPO_ROOT}/src/include ${GRACKLE_REPO_ROOT}/src/clib/autogen)
  target_link_libraries(Grackle::Grackle INTERFACE -lgrackle)
  target_link_directories(Grackle::Grackle
    INTERFACE ${GRACKLE_REPO_ROOT}/src/clib/.libs)

elseif (NOT "$ENV{PYGRACKLE_LEGACY_LINK}" STREQUAL "")

  message(FATAL_ERROR
    "when defined, the PYGRACKLE_LEGACY_LINK environment variable must be "
    "\"\" or \"classic\". It can't be: \"$ENV{PYGRACKLE_LEGACY_LINK}\"")

elseif(DEFINED ENV{Grackle_ROOT} OR DEFINED ENV{Grackle_DIR})

  # In this case, we are using CMake's find_package command to link against an
  # existing libgrackle library (built as a shared library with CMake)
  # This will work with a build-directory OR install-directory.
  # -> by default we assume that this library will NOT be moved (and will not
  #    be at a location that the system knows to search at runtime).
  # -> to accomplish this, we need to modify the install RPATH. But we only do 
  #    this if user doesn't express some other RPATH preference
  # -> ideally, would only modify the RPATH if it isn't in a standard system
  #    search path, but that is tricky to check. For now, the caller should set
  #    CMAKE_SKIP_INSTALL_RPATH=OFF in that scenario
  set(_PYGRACKLE_CONSUME_MODE "ExternalCMake")

  find_package(Grackle COMPONENTS shared REQUIRED)
  if ((NOT CMAKE_SKIP_RPATH) AND
      (NOT DEFINED CMAKE_BUILD_RPATH) AND
      (NOT DEFINED CMAKE_BUILD_RPATH_USE_ORIGIN) AND
      (NOT DEFINED CMAKE_BUILD_WITH_INSTALL_RPATH) AND
      (NOT DEFINED CMAKE_INSTALL_REMOVE_ENVIRONMENT_RPATH) AND
      (NOT DEFINED CMAKE_INSTALL_RPATH) AND
      (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH) AND
      (NOT CMAKE_SKIP_BUILD_RPATH) AND
      (NOT CMAKE_SKIP_INSTALL_RPATH))

    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH "ON")
  endif()

else()
  # standalone pygrackle-build (Grackle is freshly compiled as a shared lib
  # and packaged as part of pygrackle)
  # -> this is the most reliable/robust way to install pygrackle
  set(_PYGRACKLE_CONSUME_MODE "Embedded")

  enable_language(Fortran) # <- needed for embedded CMake builds
  set(BUILD_SHARED_LIBS ON)
  add_subdirectory(${GRACKLE_REPO_ROOT} grackle_sub-build-dir)
  # NOTE: as grackle evolves (and gets more complex), it may be better to use
  #       the built-in installation recipe. We may also want to do that based
  #       on how best to handle the conda installations.
  install(TARGETS Grackle_Grackle DESTINATION ${SKBUILD_PROJECT_NAME})
  if (APPLE)
    set(CMAKE_INSTALL_RPATH @loader_path)
  else()
    set(CMAKE_INSTALL_RPATH $ORIGIN)
  endif()
endif()

if (NOT TARGET gtest_evolve)
  add_subdirectory(${GRACKLE_REPO_ROOT}/src/grtest grtest_sub-build-dir)
endif()
# the next line may be redundant in some cases
set_target_properties(grtest_evolve PROPERTIES POSITION_INDEPENDENT_CODE True)

# Actual build Logic
# ------------------

# take the __config__.py.in template file and use it to create __config__.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/__config__.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/__config__.py @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__config__.py
  DESTINATION ${SKBUILD_PROJECT_NAME})

# transpile gracke_wrapper.pyx to a C file and then compile into a shared
# represented by the grackle_wrapper target
cython_transpile(
  pygrackle/grackle_wrapper.pyx
  LANGUAGE CXX
  CYTHON_ARGS -3
  OUTPUT_VARIABLE grackle_wrapper_cxx
)
python_add_library(grackle_wrapper MODULE "${grackle_wrapper_cxx}" WITH_SOABI)
target_compile_definitions(grackle_wrapper
  PRIVATE TRADITIONAL_IN_SOURCE_BUILD=${traditional_in_source_build})
target_link_libraries(grackle_wrapper PUBLIC Grackle::Grackle)

cython_transpile(
  pygrackle/_pygrackle_utils.pyx
  LANGUAGE CXX
  CYTHON_ARGS -3
  OUTPUT_VARIABLE _pygrackle_utils_cxx
)
python_add_library(_pygrackle_utils MODULE "${_pygrackle_utils_cxx}" WITH_SOABI)
target_link_libraries(_pygrackle_utils PUBLIC Grackle::Grackle grtest_evolve)

# Installation Logic
# ------------------
install(
  TARGETS grackle_wrapper _pygrackle_utils
  DESTINATION ${SKBUILD_PROJECT_NAME}
)


