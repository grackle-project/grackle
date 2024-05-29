# Specify Installation Rules
# ==========================
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
# -> this can be simplified after we GH PR #188 is merged
install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/clib/
  TYPE INCLUDE COMPONENT Grackle_Development
  FILES_MATCHING
  REGEX "grackle[^/]*\\.(h|def)$"
  REGEX "grackle_(macros|float).h$" EXCLUDE
  PATTERN "grackle_chemistry_data_fields.def" EXCLUDE
  REGEX "\\.libs" EXCLUDE
)

install(DIRECTORY ${GRACKLE_GENRATED_PUBLIC_HEADERS}/
  TYPE INCLUDE COMPONENT Grackle_Development
  FILES_MATCHING PATTERN "*.h" PATTERN "*.def"
)

install(TARGETS Grackle_Grackle
  # uncomment following line when ready to define GrackleTargets "export set"
  # EXPORT GrackleTargets

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


# In the short term, we have explicitly commented out the following code -- it
# deals with creating (and defining rules for installing) a "Config Package
# File".  We do this for a few reasons:
#
#  - It's not strictly necessary. The cmake build system can be used to:
#     - embed the compilation of grackle library within other cmake projects
#     - create a Grackle installation very similar to the classic libtool-based
#       build system (if you currently link against the classic installation,
#       then you can link against this installation).
#  - creating a "Config Package File" is a little complicated
#     - This is a file, that we would name GrackleConfig.cmake, that would help
#       other cmake projects link against the installed grackle libraries.
#       Among other things, it determines any compiler flags that downstream
#       applications need to link against grackle.
#     - it's a little complicated to get this right. In particular, it needs to
#       conditionally identify any transitive dependencies based on whether
#       the downstream user wants to use the shared or static library (and the
#       dependency list depends on whether grackle is compiled with openmp)
#     - it's also somewhat important to ensure that the package file provides
#       the same experience as directly embedding grackle

## generate and installs a file containing CMake code code that declares targets
## for all export-key
## - currently, only 1 target is associated with that export key
## - all targets in that file will be in the Grackle:: namespace
#install(EXPORT GrackleTargets
#  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Grackle
#  NAMESPACE Grackle::
#)
