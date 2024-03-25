# by default, a static library is built. But this can be explicitly controlled
# by externally assigning the BUILD_SHARED_LIBS variable a value of ON or OFF

# load and and set the version information
set(grackle_vfile ${CMAKE_CURRENT_SOURCE_DIR}/VERSION)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${grackle_vfile})

file(READ ${grackle_vfile} grackle_vfile_content)
if (grackle_vfile_content MATCHES
    "^([0-9]+(\\.[0-9_]+)*(\\.[-+A-Za-z0-9_]+)?)[ \t\r\n]*$")
  set(_GRACKLE_FULL_VERSION "${CMAKE_MATCH_1}")
  string(REGEX REPLACE "(\\.[-+A-Za-z0-9_]*[-+A-Za-z_][-+A-Za-z0-9_]*)" ""
    _GRACKLE_VERSION_NO_DEVTAG "${_GRACKLE_FULL_VERSION}")
else()
  message(FATAL_ERROR "Invalid version number in ${grackle_vfile}")
endif()


# Preamble
cmake_minimum_required(VERSION 3.16)
cmake_policy(SET CMP0077 NEW)
project(Grackle VERSION ${_GRACKLE_VERSION_NO_DEVTAG} LANGUAGES C Fortran)

# Specify that `./cmake` contains cmake modules defining useful functionality
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Project-Wide Configuration
include(GNUInstallDirs)

# determine whether we are building a stand-alone program, or if grackle has
# been embedded within another project
set(GRACKLE_IS_TOP_LEVEL ON)
if (NOT CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(GRACKLE_IS_TOP_LEVEL OFF)
endif()

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "Please don't perform an in-source build.")
endif()

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "no build type detected, defaulting to Release build-type")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

option(GRACKLE_USE_DOUBLE "compile the code with double precision" ON)


# Externally provided content
find_package(HDF5 REQUIRED COMPONENTS C)
add_library(GRACKLE_HDF5_C INTERFACE IMPORTED)
target_link_libraries(GRACKLE_HDF5_C INTERFACE ${HDF5_C_LIBRARIES})
target_compile_definitions(GRACKLE_HDF5_C
  INTERFACE ${HDF5_C_DEFINITIONS} -DH5_USE_16_API)
target_include_directories(GRACKLE_HDF5_C INTERFACE ${HDF5_C_INCLUDE_DIRS})

option(GRACKLE_USE_OPENMP "Use OpenMP" OFF)
message("${CMAKE_MODULE_PATH}")
if (GRACKLE_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS C Fortran)
endif()


# Main targets built by this project
add_subdirectory(src/clib)


# Packaging (not needed if there is a parent project)
if(GRACKLE_IS_TOP_LEVEL)
  # do packaging...
endif()
