# Specify Installation Rules
# ==========================

# Before we go further, we are at a cross-roads of the various scenarios we
# can assume during installation.
#
# -> These assumptions mostly come up in the context of the "linker flags"
#    a downstream application needs to specify in order to properly link against
#    each of Grackle's implicit/explicit runtime dependencies (i.e. shared
#    libraries) when using Grackle as a static library (e.g. libgrackle.a)
#
# -> We report these linker-flags to downstream applications (in one form or
#    another) in special metadata files like pkg-config's .pc files or cmake's
#    Package-Config files
#
# -> aside: when a downstream application uses Grackle as a shared library
#    (e.g. libgrackle.so), they DON'T need to specify "linker flags" for
#    Grackle's runtime dependencies
#
# For the sake of the conversation, let's imagine that the we are using `ld`
# for linking together object files and the minimal set of "linker flags" for
# Grackle's dependencies are: `-lhdf5 -lm -lgfortran`
#
# We can now consider the 3 main scenarios:
#
# 1. The simplest case: there's 1 version of each of grackle's runtime
#    dependencies (or there are no potential version incompatabilites) AND they
#    are all installed in the standard-system-locations.
#
#    - In this case, we can use `-lhdf5 -lm -lgfortran` as flags. `ld` will be
#      able to find each library at compile time, and the operating system's
#      "dynamic linker" will be able to locate each dependency at runtime
#
#    - this scenario applies to most "vanilla" unixy systems that used a
#      package manager to install the dependencies
#
# 2. In this case, we can imagine that there are multiple incompatible versions
#    of (one or more) of Grackle's runtime library dependencies and that all of
#    the versions the runtime dependencies are in standard-system-locations.
#
#    - In this case, we can use `-L` flag(s) to specify paths that encourage
#      `ld` to identify the correct versions of the runtime dependencies while
#      assembling the downstream application and will encode the correct
#      version info into the resulting binary. The operating system's "dynamic
#      linker" will then identify the correct dependency at runtime.
#
#    - We are assuming that the paths remain valid between compiling grackle
#      and compiling the downstream application.
#
#    - Alternatively, we try the `-l` flag to specify something like
#      `-l:libhdf5.so.100` to force more precise linking... (but that may not
#      be a general solution)
#
# 3. In this case, we can imagine that one or more of Grackle's runtime library
#    dependencies is found in a non-standard location.
#
#    - a solution is to specify rpaths in linker flags. We would also need to
#      modify the shared library slightly in this scenario...
#
#    - we assume that the locations of the runtime-dependencies will always be
#      unchanged whenever we try to execute the program (if anything dependency
#      ever gets moved, the downstream application will stop working). Issues
#      could also arise if a downstream application depends on a
#      backwards-compatible, newer-version of grackle's dependencies.
#
#    - there is also some ambiguity here about determining which libraries
#      to introduce rpaths for... If a dependency can currently be found by the
#      operating system's "dynamic linker", but that's only because of the
#      values currently stored in an env variable (like LD_LIBRARY_PATH), what
#      do we do? This comes up when using intel-compilers provided through
#      intel-oneAPI
#
#    - a better solution to this scenario might be to encourage people to
#      directly link against stuff in the build directory. CMake already deals
#      with these linking issues for us (in there)
#
#
# It would be nice to support all 3 case but that's complicated. The first case
# is the simplest (and most-standard case). For now, we will ONLY support that
# case. But maybe we can introduce a user-flag to support the other cases.
#
# Alternatively, the final case seems like it might be addressed if people
# simply link to stuff in the build-directory. If we reorganize the results of
# a build so that they have the same install-directory







# important note: cmake uses the term `COMPONENT` to mean different things
#
# - there are "installation components":
#   - corresponds to COMPONENT kwarg(s) appearing in cmake's `install` command
#   - these "components" are for people like a Debian packager who might use
#     them to direct cmake to "install" only subsets of files to locations that
#     would become hypothetical packages called `libgrackle` (which just
#     contains the shared library) & `libgrackle-dev` (which includes the files
#     needed to compile a program using grackle)
#   - everybody else can ignore these COMPONENTS (by default, cmake installs
#     as much as it can -- plus they have no impact on what is built)
# - there are "import components"
#   - corresponds to the `COMPONENTS` and `OPTIONAL_COMPONENTS` kwargs
#     appearing in a downstream project's `find_package` invocation
#   - these can be used to indicate that only a subset of functionality from a
#     modularized-large-project dependency are needed
#   - our plan is to do something similar to hdf5 (they use them to specify a
#     requirement on a shared or static version of a library)
#   - these `COMPONENTS` are defined within the "Config Package File" that gets
#     installed by the rule defined by invoking `install(EXPORT ...)`

# 

# define installation rules for installation of public header files
# -> this can be simplified after GH PR #188 is merged
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/include/
  TYPE INCLUDE COMPONENT Grackle_Development
  FILES_MATCHING
  REGEX ".*\\.(h|def)$" # <- needed so we don't copy any "*.h.in" files
  REGEX "grackle_float.h$" EXCLUDE # <- can be removed after GH-203 or after we
)                                  #    remove the traditional build-system

install(DIRECTORY ${GRACKLE_GENRATED_PUBLIC_HEADERS}/
  TYPE INCLUDE COMPONENT Grackle_Development
  FILES_MATCHING PATTERN "*.h" PATTERN "*.def"
)

install(TARGETS Grackle_Grackle
  # declare that the target(s) is part of the define GrackleTargets "export set"
  EXPORT GrackleTargets

  # define rules for installing all specified shared libraries (if any)
  LIBRARY
    COMPONENT Grackle_RUNTIME
    NAMELINK_SKIP
    DESTINATION ${CMAKE_INSTALL_LIBDIR} # technically, this is redundant

  # define rules for installing all specified static libraries (if any)
  ARCHIVE
    COMPONENT Grackle_Development
)

if (BUILD_SHARED_LIBS)
  # (As noted above) Because we renamed the shared library so its called
  # `libgrackle-{VERSION_NUM}.so` (rather than `libgrackle.so`), we need an
  # install-rule to make a symlink called libgrackle.so to support compilation
  # with `-lgrackle`. (This is consistent with the classic build-system)
  install(CODE "
    set(_prefix \"${CMAKE_INSTALL_PREFIX}\")
    if(DEFINED ENV{DESTDIR})
      message(WARNING
        \"linking to libgrackle.so (during install) is untested with DESTDIR\")
      set(_prefix \"\$ENV{DESTDIR}\")
    elseif(NOT IS_ABSOLUTE \${CMAKE_INSTALL_PREFIX})
      # install probably triggered by `cmake --install <p1> --prefix <p2>`
      # and <p2> is not an absolute path...
      get_filename_component(_prefix \${CMAKE_INSTALL_PREFIX} ABSOLUTE)
    endif()

    set(_COMMON \"\${_prefix}/${CMAKE_INSTALL_LIBDIR}\")
    set(_LIB \"\${_COMMON}/$<TARGET_FILE_NAME:Grackle_Grackle>\")
    set(_LINK \"\${_COMMON}/libgrackle$<TARGET_FILE_SUFFIX:Grackle_Grackle>\")

    message(STATUS \"Creating symlink to \${_LIB} called \${_LINK}\")

    execute_process(COMMAND
      ${CMAKE_COMMAND} -E create_symlink \${_LIB} \${_LINK}
    )"
  )
endif()

# precompute metadata-file information (to help with linking)
# -----------------------------------------------------------

# Interlude to define functionality for determining values that will be used
# in metadata files (e.g. a pkg-config file or cmake Package Config File)


if (GRACKLE_USE_OPENMP)
  set(_GRACKLE_OpenMP_LIBS ${OpenMP_C_LIB_NAMES})
  list(APPEND _GRACKLE_OpenMP_LIBS ${OpenMP_Fortran_LIB_NAMES})
  list(REMOVE_DUPLICATES _GRACKLE_OpenMP_LIBS)
else()
  set(_GRACKLE_OpenMP_LIBS "")
endif()

# this is used to get (non-C) language-specific implicit requirements
function(get_implicit_link_reqs lang out_libs out_linkdirs)
  set(libs ${CMAKE_${lang}_IMPLICIT_LINK_LIBRARIES})
  list(REMOVE_DUPLICATES libs)

  # remove some common dependencies (that aren't really language-specific)
  set(nondeps
    c     # libc.so     -> this is a dependency of (almost) everything
    gcc_s # libgcc_s.so -> dependency commonly needed by underlying linker
    gcc   # libgcc.a    -> dependency commonly needed by underlying linker
  )
  list(REMOVE_ITEM libs ${nondeps})
  set(${out_libs} ${libs} PARENT_SCOPE)

  set(linkdirs ${CMAKE_${lang}_IMPLICIT_LINK_DIRECTORIES})
  list(REMOVE_DUPLICATES linkdirs)
  set(${out_linkdirs} ${linkdirs} PARENT_SCOPE)
endfunction()

get_implicit_link_reqs(Fortran Fortran_implicit_libs Fortran_implicit_includes)
set(_TOOLCHAIN_LINK_LIBS ${Fortran_implicit_libs})
list(APPEND _TOOLCHAIN_LINK_LIBS m) # explicit c requirement (but may
                                    # be a duplicate)
list(REMOVE_DUPLICATES _TOOLCHAIN_LINK_LIBS)


# Define the grackle.pc file
# --------------------------

set(_GRACKLE_PC_PRIVATE_LIBS_LIST ${_TOOLCHAIN_LINK_LIBS})

if (NOT _GRACKLE_OpenMP_LIBS STREQUAL "")
  list(APPEND _GRACKLE_PC_PRIVATE_LIBS_LIST ${_GRACKLE_OpenMP_LIBS})
endif()

# we will check whether the an hdf5.pc file exists
set(_HDF5_PC_EXISTS "FALSE")
find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
  # we could be a little more careful about version numbers...
  # we explicitly avoid using the pkg-config functions provided by PkgConfig
  # -> those define a lot of variables we don't care about
  # -> we just want to check if an hdf5 pkg-config file exists in standard
  #    search paths
  execute_process(COMMAND ${PKG_CONFIG_EXECUTABLE} --modversion hdf5
    RESULT_VARIABLE _HDF5_PC_FILE_FOUND
    OUTPUT_VARIABLE _dummy_stdout_variable
    ERROR_VARIABLE _dummy_stderr_variable
  )

  if(_HDF5_PC_FILE_FOUND EQUAL 0)
    set(_HDF5_PC_EXISTS "TRUE")
  endif()
endif()

if (_HDF5_PC_EXISTS)
  set(_PC_REQUIRES_PRIVATE_LINE "Requires.private: hdf5")
else()
  set(_PC_REQUIRES_PRIVATE_LINE
    "# Requires.private omitted (hdf5.pc wasn't found)")
  list(APPEND _GRACKLE_PC_PRIVATE_LIBS_LIST "hdf5")
endif()

list(TRANSFORM _GRACKLE_PC_PRIVATE_LIBS_LIST PREPEND "-l")
string(REPLACE
  ";" " " _GRACKLE_PC_PRIVATE_LIBS "${_GRACKLE_PC_PRIVATE_LIBS_LIST}")

# retrieve list of variables conveying Grackle informational properties
include(TargetInfoProps)
get_info_properties_export_str(Grackle_Grackle
    PKG_CONFIG _GRACKLE_PC_INFO_PROPERTIES)

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/grackle.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/grackle.pc @ONLY)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/grackle.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)

# Define the cmake Package Config File
#-------------------------------------

include(CMakePackageConfigHelpers)

# The following function implements standardized logic for determining
# version compatability when a downstream project calls
#    find_project(Grackle ...)
# with a version constraint. The produced file is used to determine whether
# grackle is compatible with their requested version number (or finding an
# installed version of grackle that is compatible).
# -> Doing this "right" is a little complicated (the way dev versions were
#    named before grackle version 3.3 would have made it much harder.
# -> very modern versions of cmake (3.19 and newer) let people specify version
#    ranges. In that scenario, we could use the SameMajorVersion policy.
#    This is the policy intended for semantic versioning.
# -> earlier versions for cmake only let people supply a single version number.
#    Since we may drop deprecated functionality between minor versions, we
#    would probably want to use SameMinorVersion policy.
#
# Realistically, there won't be any issues unless you have multiple Grackle
# installations on your machine. If someone is doing that, we can assume they
# have some level of expertise. So let's just go with SameMajorVersion
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfigVersion.cmake"
  VERSION "${Grackle_VERSION}" # <- variable was set by project(Grackle ...)
  COMPATIBILITY SameMajorVersion
)

# Configure our Config-Package-File. Currently, we don't need to use
# configure_package_config_file (but we could use it in the future)

# retrieve string that we use within the template-file to define informational
# properties of the Grackle::Grackle target
get_info_properties_export_str(Grackle_Grackle
    CMAKE_CONFIG _GRACKLE_INFO_PROPERTIES)
configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/GrackleConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfig.cmake
  @ONLY
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Grackle
)

if (BUILD_SHARED_LIBS)
  set(GRACKLE_CONFIG_FILL_VAL "shared")
else()
  set(GRACKLE_CONFIG_FILL_VAL "static")
endif()

# generate and installs a file containing CMake code code that declares targets
# for all export-keys
# - currently, only 1 target is associated with that export key
# - all targets in that file will be in the Grackle:: namespace
install(EXPORT GrackleTargets
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Grackle
  NAMESPACE Grackle::
  FILE Grackle_${GRACKLE_CONFIG_FILL_VAL}_targets.cmake
)

