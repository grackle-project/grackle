cmake_minimum_required(VERSION 3.16)
project(GrackleExample LANGUAGES C)

find_package(Grackle REQUIRED)

add_executable(example_app c_example.c)

# we need to include m in the following command because c_example.c makes
# direct use of the standard C library's math functions
target_link_libraries(example_app Grackle::Grackle m)
