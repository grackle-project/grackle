# declare grtest library components
# ---------------------------------
# -> this is an internal library that defines reusable testing utilities

# we use an object library so we don't encounter problems with python bindings
add_library(grtest_evolve OBJECT
  evolve/evolve.cpp evolve/evolve.hpp
)
target_compile_features(grtest_evolve PUBLIC cxx_std_17)
target_link_libraries(grtest_evolve PUBLIC Grackle::Grackle)
# external code uses #include grtest/evolve/evolve.hpp
target_include_directories(grtest_evolve
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/../
)
