
add_executable(Grackle_Grcli
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

set_target_properties(Grackle_Grcli
  PROPERTIES CXX_STANDARD 17
  OUTPUT_NAME grcli # controls what the file is named
  #EXPORT_NAME Grcli # <- if we want to export, uncomment this
)

target_link_libraries(Grackle_Grcli PRIVATE Grackle::Grackle)

add_executable(Grackle::Grcli ALIAS Grackle_Grcli)
