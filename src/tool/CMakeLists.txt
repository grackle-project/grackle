
add_executable(grtool
  ChemistryData.h
  CliParser.h
  FieldData.h FieldData.C
  cmd_bench.h cmd_bench.C
  cmd_show.h cmd_show.C
  handle_gr_parameters.h handle_gr_parameters.C
  main.C
  utils.h utils.C
  view.h
)

set_target_properties(grtool PROPERTIES CXX_STANDARD 17)

target_link_libraries(grtool PRIVATE Grackle::Grackle)
