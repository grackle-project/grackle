
add_executable(grcli
  ChemistryData.h
  CliParser.h
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

set_target_properties(grcli PROPERTIES CXX_STANDARD 17)

target_link_libraries(grcli PRIVATE Grackle::Grackle)
