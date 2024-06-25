# This is a cmake module that defines the following
# create_example_config_Makefile function

function(create_example_config_Makefile filename target fc_flags)

  # This is a cmake function used to generate a Makefile that can be used to
  # build the various code examples.
  #
  # Parameters
  # ----------
  # filename
  #    Specifies the location where the resulting Makefile will be written
  # target
  #    Specifies the build-target corresponding to the primary shared library,
  #    that is built by the build-system. This is used by the function to
  #    determine the include and linker flags that must be used while compiling
  #    the examples
  # fc_flags
  #    These are flags that are to be passed to the fortran compiler used to
  #    compile the fortran example. This should be used to specify the general
  #    style of fortran code used throughout the codebase.
  #
  # Note
  # ----
  # This type of operation is NOT idiomatic CMake. We primarily do this for
  # supporting code-testing. This is not guarateed to work with all choices
  # of compilers.

  if (NOT BUILD_SHARED_LIBS)
    message(FATAL_ERROR
      "cannot generate Makefile for examples when building a static library "
    )
  endif()

  # when building a shared library LIBS and INCLUDES can be blank!

  file(GENERATE OUTPUT ${filename} CONTENT "
# This is an auto-generated Makefile to be used to compile the code examples
# and link against the cmake-built libgrackle shared library in the build tree
# - it's primarily intended to be used during CI
# - it may not work for all choices of compilers

# First, we specify the paths to the compilers (and to the linker) that should
# be used to build the example problems. 

CC := ${CMAKE_C_COMPILER}
CXX := $<IF:$<BOOL:${CMAKE_CXX_COMPILER}>,${CMAKE_CXX_COMPILER},/usr/bin/c++>
FC := ${CMAKE_Fortran_COMPILER}
LD := ${CMAKE_LINKER}

# Next, define any flags that should be used when compiling the examples.
# -> Ideally, these would all be blank (a shared library dependency should not
#    strongly affect how a downstream application is compiled)
# -> In practice, we use this to pass compiler-dependent fortran flags that
#    ensures that the compiler properly handles the chosen fortran dialect

CFLAGS =
CXXFLAGS =
FFLAGS := $<JOIN:${fc_flags}, >

# Here, we would identify transitive dependencies introduced by Grackle
# - for a cmake-constructed shared library, LIBS is always empty (regardless
#   of whether the example supports openmp)
# - if we constructed a static library, LIBS, would not be blank. It would
#   include linker flags for all libraries Grackle depends upon (e.g. hdf5)
# - since the grackle's public headers are all self-contained, INCLUDES is
#   always blank. If they say, included the header file from hdf5, this would
#   be non-blank

LIBS = # unnecessary for a cmake-constructed shared lib
INCLUDES = # unnecessary for a cmake-constructed shared lib

# Finally, define variables that encode information about how exactly to link
# against grackle
# -> GRACKLE_INCLUDE specifies compiler flags to specifying where to include
#    grackle's public headers from
# -> GRACKLE_LIB specifies linker flags to use for linking against the shared
#    library.
# -> special care has been taken to ensure that the flags specify that the
#    example code should link against the code currently in the build-tree
#    (rather than say a previously installed version)
# -> WARNING: these particular flags may not be portable to all compilers

GRACKLE_INCLUDE := -I$<JOIN:$<TARGET_PROPERTY:${target},INTERFACE_INCLUDE_DIRECTORIES>, -I>
GRACKLE_LIB := -Wl,-rpath,$<TARGET_FILE_DIR:${target}> -L$<TARGET_FILE_DIR:${target}> -l$<TARGET_FILE_BASE_NAME:${target}>


")

endfunction()
