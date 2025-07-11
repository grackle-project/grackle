# generate grackle_float.h AND auto_general.c
# ===========================================

# perform a check to try to help out with the scenario where the classic build
# system in an older Grackle revisions placed auto-generated files in locations
# where a cmake build might mistakenly use the wrong version of the file
# -> this is probably only an issue with grackle_float.h, but we handle all
#    auto-generated files for consistency
# -> in modern Grackle revisions, classic-build put these files in better spots
# -> I don't think this is a perfect check. It might not help if you
#      1. configure a cmake-build-directory
#      2. checkout on old grackle revision and perform a classic-build
#      3. return the latest revision and use the cmake-build-directory from 1.
#    But, this is a fairly pathological scenario
set(bad_fname_l "grackle_float.h;../include/grackle_float.h;auto_general.c")
foreach(fname IN LISTS bad_fname_l)
  if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${fname}")
    message(FATAL_ERROR
      "Please delete the ${CMAKE_CURRENT_SOURCE_DIR}/${fname} file and perform "
      "a fresh cmake-build. That file may confuse the cmake-build. Unless you "
      "manually added that file yourself, that file is left over from a "
      "classic-build in an older Grackle revision")
  endif()
endforeach()


# first, declare recipe for generating grackle_float.h:
if (GRACKLE_USE_DOUBLE)
  set(GRACKLE_FLOAT_MACRO "GRACKLE_FLOAT_8")
else()
  set(GRACKLE_FLOAT_MACRO "GRACKLE_FLOAT_4")
endif()
configure_file(../include/grackle_float.h.in
  ${GRACKLE_GENERATED_PUBLIC_HEADERS}/grackle_float.h @ONLY)

# next, declare recipe for generating auto_general.c:

# fetch necessary version information via query-version.py script
# -> GIT_BRANCH and GIT_REVISION hold sensible vals when git isn't installed
# -> VERSION_NUM may hold more information than Grackle_VERSION (e.g. like the
#    "dev" suffix)

set(_query_version
  "${CMAKE_CURRENT_SOURCE_DIR}/../../scripts/query_version.py")
function(query_version option OUTVAR)
  execute_process(
      COMMAND ${_query_version} ${option}
      RESULT_VARIABLE RSLT OUTPUT_VARIABLE OUT
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(RSLT)
    message(FATAL_ERROR "Call to `${_query_version} ${option}` failed")
  endif()
  set(${OUTVAR} ${OUT} PARENT_SCOPE)
endfunction()

query_version(show-version VERSION_NUM)
query_version(git-branch GIT_BRANCH)
query_version(git-revision GIT_REVISION)

# in the original build-system, the following string also included things like:
# include-args, linking flags, macro definitions, all other compiler flags...
string(REPLACE "\n" "\\n" SHOW_FLAGS_STR "
  CC = ${CMAKE_C_COMPILER}
  FC = ${CMAKE_Fortran_COMPILER}
  LD = ${CMAKE_LINKER}
")

string(REPLACE "\n" "\\n" SHOW_CONFIG_STR "
   Built with CMake

   GRACKLE_USE_DOUBLE                  : ${GRACKLE_USE_DOUBLE}
   GRACKLE_USE_OPENMP                  : ${GRACKLE_USE_OPENMP}
   BUILD_TYPE                          : ${CMAKE_BUILD_TYPE}
")

configure_file(auto_general.c.in auto_general.c @ONLY)

# Declare the main target
# =======================

add_library(Grackle_Grackle

  # C source files
  calculate_cooling_time.c
  calculate_dust_temperature.c
  calculate_gamma.c
  calculate_pressure.c
  calculate_temperature.c
  dynamic_api.c
  grackle_units.c
  index_helper.c
  initialize_chemistry_data.c
  initialize_cloudy_data.c
  initialize_rates.c
  initialize_UVbackground_data.c
  rate_functions.c
  set_default_chemistry_parameters.c
  solve_chemistry.c
  status_reporting.c status_reporting.h
  update_UVbackground_rates.c
  utils.c

  # auto-generated C source files
  ${CMAKE_CURRENT_BINARY_DIR}/auto_general.c

  # Fortran Source Files
  calc_tdust_1d_g.F
  calc_tdust_3d_g.F
  calc_temp1d_cloudy_g.F
  calc_temp_cloudy_g.F
  cool1d_cloudy_g.F
  cool1d_cloudy_old_tables_g.F
  cool1d_multi_g.F
  cool_multi_time_g.F
  interpolators_g.F
  solve_rate_cool_g.F

  # explicitly list the autogenerated header(s)
  # -> If we don't explicitly list these (storing it in a list isn't adequate)
  #    a subtle bug can arise in the scenario when you are freshly building
  #    Grackle some after a previous build was installed in installed in a
  #    standard system location (like /usr/local)
  # -> In that scenario, CMake doesn't realize it needs to generate the
  #    header(s) because the compiler can finds versions of the headers from
  #    the prior build in the standard system location.
  ${GRACKLE_GENERATED_PUBLIC_HEADERS}/grackle_float.h

  # although not strictly necessary, listing out other header files is a "best
  # practice" (and is needed to make them show up in certain IDEs)

  # C/C++ public headers
  ../include/grackle.h
  ../include/grackle_chemistry_data.h
  ../include/grackle_misc.h
  ../include/grackle_rate_functions.h
  ../include/grackle_types.h

  # C private headers
  cie_thin_cooling_rate_tables.h
  grackle_chemistry_data_fields.def # <-- acts as a C header
  grackle_macros.h
  index_helper.h
  phys_constants.h
  utils.h

  # Fortran public headers
  ../include/grackle.def
  ../include/grackle_fortran_interface.def
  ../include/grackle_fortran_types.def

  # Fortran private headers
  phys_const.def
)

# A downstream application built with CMake generally has 2 inclusion
# approaches (ways that they can make use of this build):
#    1. they can build Grackle as part of the downstream build.
#       - in this case, the downstream build calls add_subdirectory on the top
#         level directory of Grackle and they can access the targets defined in
#         this file.
#    2. they can use the installed library
#       - in this case, the downstream build will get access to the targets
#         defined in a special CMake export-file created during installation
#
# We take some care to make these approaches as interchangable as possible for
# downstream build. To that end, we take 2 primary steps:
#
# - we name the main target Grackle_Grackle, rather than just [Gg]rackle to
#   minimize the chance for name-collisions of targets created by the build of
#   the downstream application (this is relevant for approach #1)
#
# - Following existing conventions, the downstream project should access the
#   target called Grackle::Grackle. We make use of a namespace to follow
#   convention... (they are more useful if we built multiple targets as part of
#   the build)
add_library(Grackle::Grackle ALIAS Grackle_Grackle) # for inclusion appraoch #1

# set standard cmake-specific properties
set_target_properties(Grackle_Grackle PROPERTIES
  # specify the grackle library is properly named (when it's NOT a shared lib)
  OUTPUT_NAME "grackle"

  # choose a special name when grackle is compiled as a shared library (to
  # reduce the chance that an application linked against one version will try
  # to load a newer, incompatible version of the library)
  # -> we follow the convention of the classic build-system (and libtool) of
  #    encode encode the version in the core-part of the library name
  # -> for convenience, when Grackle is built and installed as a shared lib,
  #    custom installation logic is made to link to libgrackle.so (or
  #    libgrackle.dylib on macOS)
  LIBRARY_OUTPUT_NAME "grackle-${VERSION_NUM}"

  # specify target name that is made visible in inclusion approach #2. The
  # namespace-prefix of this target, "Grackle::", is declared in the
  # install(EXPORT ...) command included with the installation rules
  EXPORT_NAME Grackle
)

# write custom propeties used to convey information to downstream consumers
# of the Grackle library. To help us expose consistent information regardless
# of how Grackle is consume, we use this custom function that wraps
# set_target_properties (it's important to identify what is a BOOL vs STRING)
include(TargetInfoProps)
set_target_typed_info_properties(Grackle_Grackle STRING_PROPERTIES
  GRACKLE_VERSION_STR "${VERSION_NUM}"
)
set_target_typed_info_properties(Grackle_Grackle BOOL_PROPERTIES
  GRACKLE_USE_DOUBLE "${GRACKLE_USE_DOUBLE}"
  GRACKLE_USE_OPENMP "${GRACKLE_USE_OPENMP}"
)

target_include_directories(Grackle_Grackle
  # specify where to search for generated and ordinary headers when building
  # grackle AND when linking against grackle under inclusion approach #1
  # -> while it may seem unnecessary to specify the ordinary headers' directory
  #    while building grackle, it's necessary to compile auto_general.c
  PUBLIC $<BUILD_INTERFACE:${GRACKLE_GENERATED_PUBLIC_HEADERS}> # generated hdrs
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include> # public hdrs
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> # private hdrs

  # specify where to search for the other headers when linking against grackle
  # (for inclusion approach #2)
  INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(Grackle_Grackle
  PRIVATE toolchain::m
          GRACKLE_HDF5_C
          $<$<BOOL:${GRACKLE_USE_OPENMP}>:OpenMP::OpenMP_Fortran>
          $<$<BOOL:${GRACKLE_USE_OPENMP}>:OpenMP::OpenMP_C>
)

# add the necessary compiler-specific Fortran flags to ensure proper handling
# of our specific Fortran-dialect and that name mangling is performed as
# expected (these flags only get applied to Fortran source-code files)
include(AddConditionalFortranFlags)
add_conditional_fortran_flags(Grackle_Grackle)

# add the machine-specific compiler optimization options (generally specified
# in host-files) for use in the  RELEASE or RELWITHDEBINFO build-configurations
# -> as noted elsewhere, we may stop using these when Grackle isn't the
#    top-level project
# -> we explicitly add these directly to the Grackle target rather than
#    collecting them in an interface library. This lets us avoid exposing these
#    private flags when we compile Grackle_Grackle as a static library (they
#    would be exposed once we add support for the ``export`` command to allow
#    cmake projects to import Grackle from its build directory).
# -> we can start collecting them in an interface library once CMake 3.26 is our
#    minimum version, which provides the $<BUILD_LOCAL_INTERFACE:...> generator
#    expression.
target_compile_options(Grackle_Grackle PRIVATE
  "$<$<CONFIG:RELEASE>:${GRACKLE_OPTIMIZATION_FLIST}>"
  "$<$<CONFIG:RELWITHDEBINFO>:${GRACKLE_OPTIMIZATION_FLIST}>")

# define macros so that the C files know how to properly handle name-mangling
# (if cmake were our only build-system, we could use FortranCInterface module)
if ("${CMAKE_SYSTEM_NAME}" MATCHES "^(Linux)|(Darwin)$")
  target_compile_definitions(Grackle_Grackle PRIVATE "LINUX")
endif()

# define macro to suppress warnings about C files using our deprecated public
# headers (the files should include grackle.h instead). We're currently holding
# off on this to minimize conflicts with the newchem-cpp branch
target_compile_definitions(Grackle_Grackle PRIVATE GRIMPL_PUBLIC_INCLUDE=1)
