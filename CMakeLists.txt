cmake_minimum_required(VERSION 3.22)

project("GRACKLE" CXX C Fortran)

find_package(HDF5)

# Add BUILD_TESTS option available 
option(BUILD_TESTS "build tests" OFF)

# Setup GTEST if needed
if(BUILD_TESTS)
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
    include(setup_GTest) 
    add_subdirectory(tests)
endif()

add_subdirectory(src)

# Copy examples data directory
file(COPY ${CMAKE_SOURCE_DIR}/grackle_data_files DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/input DESTINATION ${CMAKE_BINARY_DIR})
