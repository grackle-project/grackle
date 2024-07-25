
add_executable(grtool
  ChemistryData.h
  main.C
  FieldData.C FieldData.h
  handle_gr_parameters.C handle_gr_parameters.h
  utils.C utils.h
  view.h
)

set_target_properties(grtool PROPERTIES CXX_STANDARD 17)

target_link_libraries(grtool PRIVATE Grackle::Grackle)
