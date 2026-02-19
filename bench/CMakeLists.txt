
add_executable(Grackle_benchmarks src/benchmarks.cpp)

set_target_properties(Grackle_benchmarks PROPERTIES
  CXX_STANDARD 17
  OUTPUT_NAME "grackle-benchmarks"  # controls what the file is named
)

target_link_libraries(Grackle_benchmarks
  PRIVATE Grackle::Grackle benchmark::benchmark grtestutils)
