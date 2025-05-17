# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(testdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)

# declare the grtest utility library
# ----------------------------------
# -> this is an internal library that defines reusable testing utilities
# -> if we add more files to this library, we should consider relocating the
#    library to a different directory

add_library(grtest_utils
  grtest_cmd.hpp grtest_cmd.cpp
  grtest_utils.hpp grtest_utils.cpp
)
# we are being a little lazy with our usage of testdeps right here
target_link_libraries(grtest_utils PUBLIC testdeps)
target_compile_features(grtest_utils PUBLIC cxx_std_17)

# these compile-definitions act as short-term hacks
target_compile_definitions(grtest_utils
  # this hack helps us get path input-file directory (we can remove it once we
  # introduce automatic file management in PR 235, PR 237, and PR 246)
  PRIVATE GR_DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/../../grackle_data_files/input/

  # this hack lets us use Operating-system specific functionality (Once PR #237
  # is merged, we should make use of the machinery introduced by that PR for
  # enabling/disabling os-specific features)
  PRIVATE "$<$<PLATFORM_ID:Linux,Darwin>:PLATFORM_GENERIC_UNIX>"
)

# start declaring targets for tests
# ---------------------------------
add_executable(runInterpolationTests test_unit_interpolators_g.cpp)
target_link_libraries(runInterpolationTests testdeps)

gtest_discover_tests(runInterpolationTests)

# one might argue that the following is more of an integration or end-to-end
# test than a unit-test
add_executable(runGhostZoneTests test_ghost_zone.cpp)
target_link_libraries(runGhostZoneTests grtest_utils testdeps)
gtest_discover_tests(runGhostZoneTests)

# this tests that the members of the chemistry_data struct can be accessed
# through the chemistry_data api (it's a special case since it involves invoking
# the system command and so it should remain separate from the other tests)
add_executable(runSyncedChemistryData test_chemistry_struct_synced.cpp)
target_link_libraries(runSyncedChemistryData grtest_utils testdeps)
target_compile_definitions(runSyncedChemistryData
  PRIVATE GR_SRC_ROOT_PATH=${PROJECT_SOURCE_DIR}
)
gtest_discover_tests(runSyncedChemistryData)
