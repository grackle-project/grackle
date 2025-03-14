# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(testdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)
target_compile_features(testdeps INTERFACE cxx_std_17)

# start declaring targets for tests
# ---------------------------------
add_executable(runFrozenKeyIdxBiMapTests test_unit_FrozenKeyIdxBiMap.cpp)
target_link_libraries(runFrozenKeyIdxBiMapTests testdeps)

gtest_discover_tests(runFrozenKeyIdxBiMapTests)


add_executable(runInterpolationTests test_unit_interpolators_g.cpp)
target_link_libraries(runInterpolationTests testdeps)

gtest_discover_tests(runInterpolationTests)


add_executable(runStatusReporting test_status_reporting.cpp)
target_link_libraries(runStatusReporting testdeps)
# short-term hack. We should use the machinery introduced in PR #237 related to
# Operating-system specific functionality when that PR is merged
target_compile_definitions(runStatusReporting PRIVATE
  "$<$<PLATFORM_ID:Linux,Darwin>:PLATFORM_GENERIC_UNIX>"
)

gtest_discover_tests(runStatusReporting)
