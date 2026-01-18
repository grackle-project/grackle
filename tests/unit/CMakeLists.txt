# declare testdeps to bundle together dependencies used by all tests
# ------------------------------------------------------------------
add_library(testdeps INTERFACE)
target_link_libraries(testdeps INTERFACE Grackle::Grackle GTest::gtest_main)

# compiler defs are short-term hacks to help tests invoke internal functions
# from Grackle (both of these will become unnecessary in the near future)
target_compile_definitions(testdeps
  # needed to invoke Fortran functions from C
  INTERFACE "$<$<PLATFORM_ID:Linux,Darwin>:LINUX>"
  # suppresses warnings about C internal headers using our deprecated public
  # headers (the files should include grackle.h instead). We're currently
  # holding off on properly fixing this to minimize conflicts with pending PRs
  # on the newchem-cpp branch
  INTERFACE GRIMPL_PUBLIC_INCLUDE=1
)

# needed to let tests call internal functions from C
target_include_directories(testdeps INTERFACE ${PROJECT_SOURCE_DIR}/src/clib)
target_compile_features(testdeps INTERFACE cxx_std_17)

# declare the grtest utility library
# ----------------------------------
# -> this is an internal library that defines reusable testing utilities
# -> if we add more files to this library, we should consider relocating the
#    library to a different directory

add_library(grtest_utils
  grtest_cmd.hpp grtest_cmd.cpp
  grtest_utils.hpp grtest_utils.cpp
  grtest_os.hpp grtest_os.cpp
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
add_executable(testFrozenKeyIdxBiMap
  test_unit_hash.cpp test_unit_FrozenKeyIdxBiMap.cpp)
target_link_libraries(testFrozenKeyIdxBiMap testdeps GTest::gmock)
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
target_link_libraries(runStatusReporting grtest_utils)
gtest_discover_tests(runStatusReporting)

# one might argue that the following is more of an integration or end-to-end
# test than a unit-test
add_executable(runVisitorTests test_visitor.cpp)
target_link_libraries(runVisitorTests testdeps)
gtest_discover_tests(runVisitorTests)

# one might argue that the following is more of an integration or end-to-end
# test than a unit-test
add_executable(runGhostZoneTests test_ghost_zone.cpp)
target_link_libraries(runGhostZoneTests grtest_utils testdeps)
gtest_discover_tests(runGhostZoneTests)

# this target tests that the members of the chemistry_data struct can be
# accessed through the "dynamic api." The test cases in this target are
# "special" since they invoke shell commands, which involves the equivalent of
# calling fork/exec or posix_spawn. There is a remote chance this could create
# problems for "death-tests", so these test-cases should remain separate from
# the rest of the gtest framework
add_executable(runSyncedChemistryData test_chemistry_struct_synced.cpp)
target_link_libraries(runSyncedChemistryData grtest_utils testdeps)
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
