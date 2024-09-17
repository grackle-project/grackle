set(GRACKLE_REPO_ROOT  ${CMAKE_CURRENT_SOURCE_DIR}/../../)

# Prologue
# --------
cmake_minimum_required(VERSION 3.16)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

# Dependencies
# ------------

# find required python components
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

# deal with dependency on Grackle and define pygrackle_build_kind variable
set(traditional_in_source_build 0) # <- will be overwritten in certain cases
if ("$ENV{PYGRACKLE_LEGACY_LINK}" STREQUAL "classic")

  # match old setuptools-approach for the builds with the classic build-system
  # -> explicitly match legacy behavior: link against library in local build &
  #    HOPE that the library will be installed when we need it
  # -> more idiomatic alternative: provide cmake with full library path
  set(traditional_in_source_build 1)
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
  set(traditional_in_source_build 0 PARENT_SCOPE)
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
  # -> this is the most reliable way to do it get pygrackle
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

# Actual build Logic
# ------------------
# -> we build a shared-library represented by the grackle_wrapper library

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/__config__.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/__config__.py @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/__config__.py
  DESTINATION ${SKBUILD_PROJECT_NAME})
  

add_custom_command(
  OUTPUT grackle_wrapper.c
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/grackle_wrapper.c from ${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/grackle_wrapper.pyx"
  COMMAND Python::Interpreter -m cython
          -3
          "${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/grackle_wrapper.pyx"
          --output-file grackle_wrapper.c
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/grackle_wrapper.pyx
          ${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/grackle_defs.pxd
  VERBATIM)

python_add_library(grackle_wrapper MODULE grackle_wrapper.c WITH_SOABI)
target_compile_definitions(grackle_wrapper
  PRIVATE TRADITIONAL_IN_SOURCE_BUILD=${traditional_in_source_build})
target_link_libraries(grackle_wrapper PUBLIC Grackle::Grackle)
target_include_directories(grackle_wrapper PRIVATE ${Python_NumPy_INCLUDE_DIRS})

# Installation Logic
# ------------------
install(TARGETS grackle_wrapper DESTINATION ${SKBUILD_PROJECT_NAME})

