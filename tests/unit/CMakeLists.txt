# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(testdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)
target_compile_features(testdeps INTERFACE cxx_std_17)

# short-term hack to get path input-file directory
# - an alternative, slightly nicer hack is introduced by PR #254
# - we can totally remove this hack once we introduce automatic file management
#   in PR 235, PR 237, and PR 246
target_compile_definitions(testdeps
  INTERFACE GR_DATADIR=${CMAKE_CURRENT_SOURCE_DIR}/../../grackle_data_files/input/
)

# start declaring targets for tests
# ---------------------------------
add_executable(runInterpolationTests
  test_unit_interpolators_g.cpp test_interpolators_comparisons.cpp
)
target_link_libraries(runInterpolationTests testdeps)

gtest_discover_tests(runInterpolationTests)
