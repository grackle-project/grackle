add_executable(runLocalSolveChemistryTest test_local_solve_chemistry.cpp)
target_link_libraries(runLocalSolveChemistryTest grackle gtest gtest_main hdf5)
target_include_directories(runLocalSolveChemistryTest PRIVATE ${PROJECT_SOURCE_DIR}/src/clib ${PROJECT_SOURCE_DIR}/src/clib/generated)
add_test(NAME LocalSolveChemistryMain COMMAND runLocalSolveChemistryTest)


