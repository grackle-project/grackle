set(GRACKLE_ROOT_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/../../)

# here we define a function that we use for backwards compatability...
# hopefully we can remove it in the future!
function(legacy_find_grackle)

  add_library(Grackle::Grackle INTERFACE IMPORTED)
  target_include_directories(Grackle::Grackle
    INTERFACE ${GRACKLE_ROOT_DIR}/src/clib ${GRACKLE_ROOT_DIR}/src/include)

  # it would be more idiomatic to provide cmake with the full path of the library.
  # but, we explicitly match legacy behavior of linking against the library in the
  # local build and HOPING that the library will be installed when we need it
  # (thats why we target_link_directories)
  target_link_libraries(Grackle::Grackle INTERFACE -lgrackle)

  if($ENV{PYGRACKLE_LEGACY_LINK} STREQUAL "classic") # traditional build
    set(traditional_in_source_build 1 PARENT_SCOPE)
    target_link_directories(Grackle::Grackle
      INTERFACE ${GRACKLE_ROOT_DIR}/src/clib/.libs)

  elseif($ENV{PYGRACKLE_LEGACY_LINK} MATCHES "cmake:(.+)")
    # we can remove this branch as soon as we figure out how to better deal with
    # the test_code_examples unit test
    set(traditional_in_source_build 0 PARENT_SCOPE)
    target_link_directories(Grackle::Grackle INTERFACE "${CMAKE_MATCH_1}/src/clib")
    target_include_directories(Grackle::Grackle
      INTERFACE "${CMAKE_MATCH_1}/generated_public_headers")
  else()
    message(FATAL_ERROR
      "when defined, the PYGRACKLE_LEGACY_LINK environment variable must be "
      "\"\", \"classic\", or  \"cmake:path/to/build\"")
  endif()
endfunction()

# Prologue
# --------
cmake_minimum_required(VERSION 3.16)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C)

# Dependencies
# ------------

# find required python components (we currently need NumPy's include directory)
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# deal with dependency on Grackle and define pygrackle_build_kind variable
set(traditional_in_source_build 0) # <- will be overwritten in certain cases
if (DEFINED ENV{PYGRACKLE_LEGACY_LINK} AND 
    NOT (ENV{PYGRACKLE_LEGACY_LINK} STREQUAL ""))
  legacy_find_grackle()
else()
  # standalone pygrackle-build (Grackle is freshly compiled as part of pygrackle)
  # -> at the moment, we compile Grackle as a static library, with position
  #    independent code (so that we only produce 1 shared library at the end)
  # -> longer term, it might be more sensible to compile it as a shared library and
  #    package that side-by-side with the python-module shared lib (that way, we
  #    could potentially provide it in conda-installations)
  enable_language(Fortran)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # 2nd arg gives build-subdirectory, where intermediate build products are put
  add_subdirectory(${GRACKLE_ROOT_DIR} grackle_sub-build)
endif()

# Actual build Logic
# ------------------
# -> we build a shared-library represented by the grackle_wrapper library

add_custom_command(
  OUTPUT grackle_wrapper.c
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/grackle_wrapper.c from ${CMAKE_CURRENT_SOURCE_DIR}/pygrackle/grackle_wrapper.pyx"
  COMMAND Python::Interpreter -m cython
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

