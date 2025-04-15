
# overwrite the global variable that controls where the binaries compiled in
# this directory are placed
# -> because this file is only read by CMake via add_subdirectory and we are
#    are not making a CACHE variable, this won't affect anything outside of
#    this directory
# -> for a build directory, <build>, this places the binaries at 
#    <build>/examples
# -> assuming that the build-directory is at the root level of the grackle
#    repository, then executing the examples from within the examples directory
#    should ensure that the paths to the datafiles are correct
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../examples)

if (GRACKLE_USE_DOUBLE)
  # this example requires C99 or newer since it uses designated initializers
  add_executable(c99_simple_example c99_simple_example.c)
  target_link_libraries(c99_simple_example Grackle::Grackle)
  target_compile_features(c99_simple_example PUBLIC c_std_99)
endif()

add_executable(c_example c_example.c)
target_link_libraries(c_example Grackle::Grackle)

add_executable(c_local_example c_local_example.c)
target_link_libraries(c_local_example Grackle::Grackle)

set(CMAKE_CXX_STANDARD 11) # needed by apple-clang compiler for some examples

add_executable(cxx_example cxx_example.C)
target_link_libraries(cxx_example Grackle::Grackle)

if (GRACKLE_USE_OPENMP)
  add_executable(cxx_omp_example cxx_omp_example.C)
  target_link_libraries(cxx_omp_example Grackle::Grackle)
endif()

# the following target uses -O1 optimization because problems can arise with
# higher levels of optimization. A comment in fortran_example.F provides more
# details about the problems
add_executable(fortran_example fortran_example.F)
target_compile_options(fortran_example PRIVATE "-O1")
target_link_libraries(fortran_example Grackle::Grackle)
