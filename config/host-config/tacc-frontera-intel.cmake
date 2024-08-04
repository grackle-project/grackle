#=======================================================================
# Provides host-file settings for frontera using the intel compilers
#
# About host-files:
# - a host-file provides machine-specific configuration options to help compile
#   Grackle as a standalone library on common HPC platforms. They should not
#   generally be used when Grackle is embedded within the build of another
#   cmake project.
#
# - These mostly exist for convenience and parity with the machine-files used
#   by the classic build-system. They are most useful on HPC systems that
#   provide multiple compiler toolchains. (They usually aren't needed on local
#   systems)
#
# - In terms of modern, idiomatic CMake features, a host-file could be
#   replaced by a combination of a toolchain-file and a preset-file
#   -> toolchain files usually define compiler-toolchain related-options and
#      are commonly used for cross-compiling. Rule of thumb: you should be able
#      to recycle toolchain-files between unrelated projects (i.e. they don't
#      include project-specific variables).
#   -> a preset file (``CMakePresets.json`` and ``CMakeUserPresets.json``) are
#      intended to be used to specify common project-specific compilation
#      options. These can be read by IDEs.
#   -> after we update the minimum required CMake version for compiling Grackle
#      to at least 3.19, we may transition to using these features.
#
# All configuration options specified in this file should use CMake's ``set``
# command with the ``CACHE`` option.
#=======================================================================


#-----------------------------------------------------------------------
# Specify Compilers:
#-----------------------------------------------------------------------
set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "")


#-----------------------------------------------------------------------
# Optimization Compiler Options:
# -> these should really be machine-specific
# -> these flags won't be used in debug-builds
# -> See the installation documentation for more details on this topic
#    (including how to properly quote groups of variables)
#-----------------------------------------------------------------------
# The Frontera User Guide recommends the following flag for their Cascade Lake
# (CLX) Compute Nodes
set(GRACKLE_OPTIMIZATION_FLIST_INIT "-xCORE-AVX512" CACHE STRING "")

#-----------------------------------------------------------------------
# Other Options:
# -> on certain platforms, if HDF5 is found in an extremely atypical spot, you
#    you might place that hint here
# -> If you find that you NEED to configure other options in order to get
#    Grackle to successfully compile (e.g. a fortran flag to tell the compiler
#    how to handle Grackle's Fortran dialect), please open an issue on GitHub.
#    These sort of issues should be handled by logic executed in Grackle's
#    CMakeLists.txt files so that we don't duplicate this kind of logic across
#    host-files of different machines (and so that Grackle can be embedded 
#    within other CMake projects)
#-----------------------------------------------------------------------

# the following frontera-specific snippet is only included for the sake of
# example (it definitely isn't required)
#   Here we provide a hint about where to find hdf5 (cmake is usually smart
#   enough that it can find it without the hint). When you execute
#       module spider hdf5/1.x.y
#   The output informs us that loading that lmod-module will define the
#   TACC_HDF5_DIR environment variable, which specifies the location
#   of the module's hdf5 installation. Here, we're telling CMake that it MUST
#   use this particular installation
set(HDF5_ROOT "$ENV{TACC_HDF5_DIR}" CACHE STRING
  "HDF5 install-location based on the loaded hdf5-module during configuration")
