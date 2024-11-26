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

add_executable(runUnitCalcTDust1DTests test_unit_calc_tdust_1d_g.cpp)
target_link_libraries(runUnitCalcTDust1DTests grackle gtest gtest_main hdf5)

gtest_discover_tests(runUnitCalcTDust1DTests)

add_executable(runUnitCool1DCloudyTests test_unit_cool1d_cloudy_g.cpp)
target_link_libraries(runUnitCool1DCloudyTests grackle gtest gtest_main hdf5)

gtest_discover_tests(runUnitCool1DCloudyTests)

add_executable(runUnitCalcTempUtilitiesTest test_unit_calc_temp_utilities_g.cpp)
target_link_libraries(runUnitCalcTempUtilitiesTest grackle gtest gtest_main hdf5)

gtest_discover_tests(runUnitCalcTempUtilitiesTest)

