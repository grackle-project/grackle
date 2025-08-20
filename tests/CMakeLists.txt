add_subdirectory(unit)

# down below, we add tests for various code-examples

# determine the path to the script used to assist with tests of code examples
set(checker_script ${CMAKE_CURRENT_LIST_DIR}/scripts/code_example_checker.py)

# Dealing with hard-coded paths to data-files
# -------------------------------------------
# These tests are complicated by the fact that code-examples all hard-code
# (relative) paths to data-files.
#
# Here we implement a 2-part stop-gap solution to gracefully handle these cases
# (The rigorously correct solution will be possible once PRs #235, #237, and
# #246 are merged into Grackle)
#
# Part 1) we launch all of the tests from the source-directory including the
# source-code files
# -> while we tried to configure the build system in such a way that the
#    code-example binaries can be executed from the directory where the
#    binaries are placed. But, this makes 2 assumptions:
#    (1) that Grackle is the top-level project
#    (2) the Grackle build directory's parent directory is the root-directory
#        of the Grackle repository
# -> executing the code-examples from the src-directory during testing will be
#    more rigorous. It should work when Grackle is compiled as part of a larger
#    project
set(example_exec_dir ${PROJECT_SOURCE_DIR}/src/example)
#
# Part 2) we perform a prerequisite test confirming that the assumed data file
# actually exists (i.e. if Grackle's submodules got cloned).
# -> we use ctest's fixture machinery to gracefully skip other tests if this
#    one fails
# -> as a bonus, this confirms that the script used to run most of these tests
#    is working properly

add_test(
  NAME GrExample.SatisfyPrereqs
  COMMAND ${checker_script} datacheck
    --exec-dir ${example_exec_dir}
    --assumed-data-path ../../input/CloudyData_UVB=HM2012.h5
)
set_tests_properties(GrExample.SatisfyPrereqs PROPERTIES
  FIXTURES_SETUP GREXAMPLE_PREREQS
)

# this case simply runs the example with 4 openmp threads
if (GRACKLE_USE_OPENMP)
add_test(NAME GrExample.cmp_omp
  COMMAND cxx_omp_example
  WORKING_DIRECTORY ${example_exec_dir}
)
set_tests_properties(GrExample.cmp_omp PROPERTIES
  ENVIRONMENT "OMP_NUM_THREADS=4"
  FIXTURES_REQUIRED GREXAMPLE_PREREQS
  DISABLED ${GRACKLE_USE_OPENMP}
)
endif()

# the remaining tests are more rigorous:
#
# first, we check the results of c_example against the historical results
# (which we previously recorded to a JSON file)
# -> we may need to loosen up rtol so that this works across machines
# -> if we ever need to regenerate the reference JSON file, we can run
#    `code_example_checker.py genjson --target <path> --output <path>`
add_test(
  NAME GrExampleHistoricalCMP.c
  COMMAND ${checker_script} cmp
    --target $<TARGET_FILE:c_example>
    --ref ${CMAKE_CURRENT_LIST_DIR}/scripts/c_example_ref.json
    --rtol 2e-12
    --exec-timeout 10
    --exec-dir ${example_exec_dir}
)
set_tests_properties(GrExampleHistoricalCMP.c PROPERTIES
  FIXTURES_REQUIRED GREXAMPLE_PREREQS
)

# the remaining tests compare the results printed by various code examples
# against each other. In general, there should be no differences in the
# results
add_test(
  NAME GrExampleCMP.c_local
  COMMAND ${checker_script} cmp
    --target $<TARGET_FILE:c_local_example>
    --ref $<TARGET_FILE:c_example>
    --exec-timeout 10
    --exec-dir ${example_exec_dir}
)
set_tests_properties(GrExampleCMP.c_local PROPERTIES
  FIXTURES_REQUIRED GREXAMPLE_PREREQS
)

add_test(
  NAME GrExampleCMP.cxx
  COMMAND ${checker_script} cmp
    --target $<TARGET_FILE:cxx_example>
    --ref $<TARGET_FILE:c_example>
    --exec-timeout 10
    --exec-dir ${example_exec_dir}
)
set_tests_properties(GrExampleCMP.cxx PROPERTIES
  FIXTURES_REQUIRED GREXAMPLE_PREREQS
)

if (TARGET fortran_example)
  # we allow for the Fortran to deviate from c_example to a small degree
  # -> to require an exact match, we probably need to rewrite all of the
  #    examples so that they use hardcoded values
  add_test(
    NAME GrExampleCMP.fortran
    COMMAND ${checker_script} cmp
      --target $<TARGET_FILE:fortran_example>
      --ref $<TARGET_FILE:c_example>
      --rtol 6e-7
      --exec-timeout 10
      --exec-dir ${example_exec_dir}
  )
  set_tests_properties(GrExampleCMP.fortran PROPERTIES
    FIXTURES_REQUIRED GREXAMPLE_PREREQS
  )
endif()

