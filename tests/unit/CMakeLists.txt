# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(testdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)

# start declaring targets for tests
# ---------------------------------
add_executable(runInterpolationTests
  test_unit_interpolators_g.cpp interpolators_g.F
)
target_link_libraries(runInterpolationTests testdeps)

add_executable(runInterpolationTests2
  test_interpolators.C utest_helpers.hpp interpolators_g.F
)
target_link_libraries(runInterpolationTests2 testdeps)

gtest_discover_tests(runInterpolationTests)
