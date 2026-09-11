cmake_minimum_required(VERSION 3.16)
# we need to enable C, to make a call to `find_package(hdf5)` work (this is
# called within `find_package(Grackle REQUIRED)`
project(GrackleExample LANGUAGES C Fortran)

find_package(Grackle REQUIRED)

add_executable(example_app fortran_example.F)

target_link_libraries(example_app Grackle::Grackle)
