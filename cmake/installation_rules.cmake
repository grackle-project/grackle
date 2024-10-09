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

include("InstallationHelperFuncs")

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

include(CreateProgram-grdata)
# define installation rules for installation of grdata cli tool. This is
# analogous to install(TARGETS ...), but is needed since grdata isn't a
# compiled executable
# -> in the future, maybe we should make this part of the Grackle_Development
#    "installation component?" (from my perspective it's probably better err
#    on the side of being a little too atomic here)
grdata_install_export_helper_(INSTALL_TOOL
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT Grackle_Tools
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

get_implicit_link_reqs(Fortran Fortran_implicit_libs Fortran_implicit_linkdirs)
set(_TOOLCHAIN_LINK_LIBS ${Fortran_implicit_libs})

# on most unix-like platforms (but not macOS), we need to explicitly link to
# to the standard library's math functions
# -> here we determine based on whether our custom toolchain::m target
#    is a dummy placeholder or not whether to add this target
get_target_property(toolchain_m_prop toolchain::m IMPORTED_LIBNAME)
if(${toolchain_m_prop})
  list(APPEND _TOOLCHAIN_LINK_LIBS m) # explicit c requirement (but may
                                      # be a duplicate)
endif()

list(REMOVE_DUPLICATES _TOOLCHAIN_LINK_LIBS)


# Define the grackle.pc file
# --------------------------

set(_STATIC_EXTRA_LINK_LIBS ${_TOOLCHAIN_LINK_LIBS})

if (NOT _GRACKLE_OpenMP_LIBS STREQUAL "")
  list(APPEND _STATIC_EXTRA_LINK_LIBS ${_GRACKLE_OpenMP_LIBS})
endif()

# we will check whether the hdf5.pc file exists
check_h5pc_exists("${HDF5_VERSION}" _HDF5_PC_EXISTS)

if (_HDF5_PC_EXISTS AND ("${HDF5_VERSION}" STREQUAL ""))
  set(_PC_STATIC_REQUIRES "hdf5")
elseif(_HDF5_PC_EXISTS)
  set(_PC_STATIC_REQUIRES "hdf5 = ${HDF5_VERSION}")
else()
  list(APPEND _STATIC_EXTRA_LINK_LIBS "hdf5")
endif()

list(TRANSFORM _STATIC_EXTRA_LINK_LIBS PREPEND "-l")
string(REPLACE
  ";" " " _STATIC_EXTRA_LINK_LIBS "${_STATIC_EXTRA_LINK_LIBS}")

# at the moment, the only extra library search paths that we are specifying is
# for finding implicit Fortran dependencies
# -> this probably isn't adequate on systems where hdf5.pc can't be found
# -> it may not be adequate on some systems when using OpenMP
#
# Our current solution currently has some warts:
# -> On most systems these flags only tell the linker where to search for the 
#    libraries at link-time (they generally don't impact the runtime search).
#    They also say to use the libraries shipped with a compiler at link time
#    -> we effectively assume that libraries are found at run-time in normal
#       system installation paths. Thus, if a downstream is linked against a
#       static grackle library it will hopefully keep working even if we
#       replace our system's Fortran compiler
#    -> it kinda makes sense to use libraries shipped with the Fortran compiler
#       for linking. If there are multiple versions of the Fortran runtime, this
#       ensures that the downstream application knows to use the runtime that
#       is compatible with the original compiler (not sure if this is really
#       an issue in-practice...). At the same time, our pkg-config file will
#       break if we ever replace/upgrade the compiler...
#    -> a compromise: link against an ABI compatible version of the runtime. If
#       we know the library exists in a standard search-path at compile-time,
#       use that instead...
# -> We also aren't very consistent. It turns out on macOS, the linker uses the
#    link-time locations at runtime as well. This is equivalent to us also
#    specifying -rpath with absolute-paths on most systems.
#    -> it turns out that this works to our advantage right now because
#       libgfortran isn't at a system install-path... it is only attached to
#       the version of libgfortran shipped with the compiler
#    -> this does mean that an installation could break if you remove/replace
#       your fortran compiler. We could potentially reduce the chance of
#       breakage during gfortran upgrades by replacing the -L path to make use
#       of the symlink at /opt/homebrew/lib/gcc/...


set(_STATIC_EXTRA_LINK_DIRS ${Fortran_implicit_linkdirs})
list(TRANSFORM _STATIC_EXTRA_LINK_DIRS PREPEND "-L")
string(REPLACE
  ";" " " _STATIC_EXTRA_LINK_DIRS "${_STATIC_EXTRA_LINK_DIRS}")

set(_GRACKLE_PC_PRIVATE_LIBS
  "${_STATIC_EXTRA_LINK_DIRS} ${_STATIC_EXTRA_LINK_LIBS}")

# retrieve list of variables conveying Grackle informational properties
include(TargetInfoProps)
get_info_properties_export_str(Grackle_Grackle
  PKG_CONFIG _GRACKLE_PC_INFO_PROPERTIES)


set(INSTALL_METADATA_DIR "${CMAKE_CURRENT_BINARY_DIR}/install-metadata") 

foreach(suffix IN ITEMS "conventional.pc" "static.pc")
  set(_extra_arg "")
  if(${suffix} STREQUAL "static.pc")
    set(_extra_arg "STATIC_ONLY")
  endif()

  configure_pkgconfig_file(
    DESTINATION ${INSTALL_METADATA_DIR}/grackle-${suffix}
    STATIC_LIBS "${_STATIC_EXTRA_LINK_DIRS} ${_STATIC_EXTRA_LINK_LIBS}"
    STATIC_REQUIRES "${_PC_STATIC_REQUIRES}"
    INFO_PROPERTIES "${_GRACKLE_PC_INFO_PROPERTIES}"
    ${_extra_arg})
endforeach()

if (BUILD_SHARED_LIBS)
  install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/grackle-conventional.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    RENAME grackle.pc)
else()
  # if shared library was previously installed, install grackle-conventional.pc
  # as grackle.pc. Otherwise, install grackle-static.pc as grackle.pc
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
    set(_PCDIR \"\${_COMMON}/pkgconfig\")
    set(_DEST \"\${_PCDIR}/grackle.pc\")
    if ((EXISTS \"\${_COMMON}/libgrackle.so\") OR
        (EXISTS \"\${_COMMON}/libgrackle.dylib\"))
      set(_SRC \"${INSTALL_METADATA_DIR}/grackle-conventional.pc\")
      file(REMOVE \${_DEST})
    else()
      set(_SRC \"${INSTALL_METADATA_DIR}/grackle-static.pc\")
    endif()

    message(STATUS \"Copying \${_SRC} to \${_DEST}\")

    execute_process(COMMAND ${CMAKE_COMMAND} -E copy \${_SRC} \${_DEST})"
  )
endif()

# Define the cmake Package Config File
#-------------------------------------

# create variable that holds the path to the current directory
set(LOCAL_CMAKE_MODULE_DIR "${CMAKE_CURRENT_LIST_DIR}")

# create variable storing where copies of the cmake files are stored so that
# they can be used without a full installation
set(BUILDTREE_CMAKE_DIR ${GRACKLE_BUILD_EXPORT_PREFIX_PATH}/cmake/Grackle)

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
  ${LOCAL_CMAKE_MODULE_DIR}/GrackleConfig.cmake.in
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

# call the analog of install(EXPORT ...) for the grdata tool
grdata_install_export_helper_(INSTALL_EXPORT
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Grackle
  FILE Grackle_grdata_pseudotarget.cmake
  TMPDIR ${CMAKE_CURRENT_BINARY_DIR}/install-metadata
)


# copy some files to the build tree to locations specified by the
# BUILDTREE_CMAKE_DIR variable so so that external cmake projects can use
# find_package to directly import Grackle::Grackle from the build-tree (without
# requiring a full installation)


file(COPY
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/install-metadata/GrackleConfigVersion.cmake
  DESTINATION ${BUILDTREE_CMAKE_DIR}
)

export(EXPORT GrackleTargets
  FILE ${BUILDTREE_CMAKE_DIR}/Grackle_${GRACKLE_CONFIG_FILL_VAL}_targets.cmake
  NAMESPACE Grackle::
)

# analog to export(EXPORT ...) for the grdata tool
grdata_install_export_helper_(EXPORT_EXPORT
  FILE ${BUILDTREE_CMAKE_DIR}/Grackle_grdata_pseudotarget.cmake
)
