
include(FetchContent)
FetchContent_Declare(
  benchmark
  # the following url contains the hash for version 1.8.5
  URL https://github.com/google/benchmark/archive/a6ad7fbbdc2e14fab82bb8a6d27760d700198cbf.zip
)

set (BENCHMARK_ENABLE_GTEST_TESTS OFF)
FetchContent_MakeAvailable(benchmark)


add_executable(Grackle_Grcli
  ChemistryData.h
  args/CliParser.h
  FieldData.h FieldData.C
  cmd_bench.h cmd_bench.C
  cmd_show.h cmd_show.C
  handle_gr_parameters.h handle_gr_parameters.C
  grid_problem.h grid_problem.C
  main.C
  operation.h operation.C
  utils.h utils.C
  view.h
)

set_target_properties(Grackle_Grcli
  PROPERTIES CXX_STANDARD 17
  OUTPUT_NAME grcli # controls what the file is named
  #EXPORT_NAME Grcli # <- if we want to export, uncomment this
)

target_include_directories(Grackle_Grcli
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(Grackle_Grcli
  PRIVATE Grackle::Grackle benchmark::benchmark)

target_compile_definitions(Grackle_Grcli PUBLIC USE_BENCHMARK=1)

add_executable(Grackle::Grcli ALIAS Grackle_Grcli)
