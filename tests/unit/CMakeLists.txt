# NOTE: the testdeps target is defined within ../grtestutils/CMakeLists.txt

add_executable(testFrozenKeyIdxBiMap test_unit_FrozenKeyIdxBiMap.cpp)
target_link_libraries(testFrozenKeyIdxBiMap testdeps)
gtest_discover_tests(testFrozenKeyIdxBiMap)

add_executable(runInterpolationTests test_unit_interpolators_g.cpp)
target_link_libraries(runInterpolationTests testdeps)
gtest_discover_tests(runInterpolationTests)

add_executable(runLinAlgTests test_linalg.cpp)
target_link_libraries(runLinAlgTests testdeps)
gtest_discover_tests(runLinAlgTests)

add_executable(runGrainSpeciesIntoTests test_grain_species_info.cpp)
target_link_libraries(runGrainSpeciesIntoTests testdeps)
gtest_discover_tests(runGrainSpeciesIntoTests)

add_executable(runStatusReporting test_status_reporting.cpp)
target_link_libraries(runStatusReporting testdeps)
gtest_discover_tests(runStatusReporting)

# one might argue that the following is more of an integration or end-to-end
# test than a unit-test
add_executable(runVisitorTests test_visitor.cpp)
target_link_libraries(runVisitorTests testdeps)
gtest_discover_tests(runVisitorTests)

# tests of the API functions
# -> one might argue that these are better classified as integration or
#    end-to-end tests than as unit-test
add_executable(runApiTests test_ghost_zone.cpp test_api_ratequery.cpp)
target_link_libraries(runApiTests testdeps)
gtest_discover_tests(runApiTests)

# this target tests that the members of the chemistry_data struct can be
# accessed through the "dynamic api." The test cases in this target are
# "special" since they invoke shell commands, which involves the equivalent of
# calling fork/exec or posix_spawn. There is a remote chance this could create
# problems for "death-tests", so these test-cases should remain separate from
# the rest of the gtest framework
add_executable(runSyncedChemistryData test_chemistry_struct_synced.cpp)
target_link_libraries(runSyncedChemistryData testdeps)
target_compile_definitions(runSyncedChemistryData
  PRIVATE
    READER_PATH=${PROJECT_SOURCE_DIR}/tests/scripts/castxml_output_reader.py
    XML_PATH=${CMAKE_CURRENT_BINARY_DIR}/grackle.h.xml
)
# every (re)build runSyncedChemistryData triggers this custom command:
# -> the castxml_wrapper.py script forwards all arguments to the castxml tool,
#    if the tool is installed. Otherwise, it injects an empty file at the output
#    location
# -> the information in the XML-file is always produced from the same version of
#    the header file used to compile runSyncedChemistryData (even if the header
#    file is modified)
add_custom_command(
  TARGET runSyncedChemistryData POST_BUILD
  COMMAND ${PROJECT_SOURCE_DIR}/tests/scripts/castxml_wrapper.py
                  -c
                  -I${GRACKLE_GENERATED_PUBLIC_HEADERS}
                  -x c++
                  --castxml-cc-gnu g++
                  --castxml-output=1
                  -o "grackle.h.xml"
                  ${PROJECT_SOURCE_DIR}/src/include/grackle.h
  BYPRODUCTS "grackle.h.xml"
  VERBATIM
)

gtest_discover_tests(runSyncedChemistryData)
