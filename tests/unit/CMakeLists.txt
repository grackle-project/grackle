# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(unittestdeps INTERFACE)
target_link_libraries(unittestdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(unittestdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(unittestdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)

# start declaring targets for tests
# ---------------------------------
add_executable(runInterpolationTests test_unit_interpolators_g.cpp)
target_link_libraries(runInterpolationTests unittestdeps)

gtest_discover_tests(runInterpolationTests)
