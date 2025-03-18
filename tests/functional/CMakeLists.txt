# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(functestdeps INTERFACE)
target_link_libraries(functestdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# short-term hack to let tests invoke Fortran functions from C
target_compile_definitions(functestdeps INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>")
target_include_directories(functestdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)

# start declaring targets for tests
# ---------------------------------
add_executable(runLocalSolveChemistryTest test_local_solve_chemistry.cpp)
target_link_libraries(runLocalSolveChemistryTest functestdeps)

gtest_discover_tests(runLocalSolveChemistryTest)
