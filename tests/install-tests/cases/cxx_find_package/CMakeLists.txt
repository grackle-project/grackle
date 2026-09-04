cmake_minimum_required(VERSION 3.16)
project(GrackleExample LANGUAGES CXX)

find_package(Grackle REQUIRED)

add_executable(example_app cxx_example.cpp)
target_link_libraries(example_app Grackle::Grackle)
