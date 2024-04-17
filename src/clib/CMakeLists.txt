set(C_SOURCES
  calculate_cooling_time.c
  calculate_dust_temperature.c
  calculate_gamma.c
  calculate_pressure.c
  calculate_temperature.c
  dynamic_api.c
  grackle_units.c
  index_helper.c
  initialize_chemistry_data.c
  initialize_cloudy_data.c
  initialize_rates.c
  initialize_UVbackground_data.c
  rate_functions.c
  set_default_chemistry_parameters.c
  solve_chemistry.c
  update_UVbackground_rates.c
  calculate_cooling_time.c
  calculate_dust_temperature.c
  calculate_gamma.c
  calculate_pressure.c
  calculate_temperature.c
  dynamic_api.c
  grackle_units.c
  index_helper.c
  initialize_chemistry_data.c
  initialize_cloudy_data.c
  initialize_rates.c
  initialize_UVbackground_data.c
  rate_functions.c
  set_default_chemistry_parameters.c
  solve_chemistry.c
  update_UVbackground_rates.c
  calculate_cooling_time.c
  calculate_dust_temperature.c
  calculate_gamma.c
  calculate_pressure.c
  calculate_temperature.c
  dynamic_api.c
  grackle_units.c
  index_helper.c
  initialize_chemistry_data.c
  initialize_cloudy_data.c
  initialize_rates.c
  initialize_UVbackground_data.c
  rate_functions.c
  set_default_chemistry_parameters.c
  solve_chemistry.c
  update_UVbackground_rates.c
  utils.c)

set(C_WRAPPER_SOURCES
  c_wrappers/wrap_interpolators_g.c)

set(C_GENERATED_SOURCES
  generated/auto_show_config.c
  generated/auto_show_flags.c
  generated/auto_show_version.c)

set(F_SOURCES
  calc_tdust_1d_g.F
  calc_tdust_3d_g.F
  calc_temp1d_cloudy_g.F
  calc_temp_cloudy_g.F
  cool1d_cloudy_g.F
  cool1d_cloudy_old_tables_g.F
  cool1d_multi_g.F
  cool_multi_time_g.F
  interpolators_g.F
  solve_rate_cool_g.F)

#set(C_DEPS_SOURCES
#  ${PROJECT_SOURCE_DIR}/src/deps/auto_show_config.c
#  ${PROJECT_SOURCE_DIR}/src/deps/auto_show_flags.c
#  ${PROJECT_SOURCE_DIR}/src/deps/auto_show_version.c)

set(C_INCLUDES_DIRECTORY ${PROJECT_SOURCE_DIR}/src/clib ${PROJECT_SOURCE_DIR}/src/clib/generated)

add_library(grackle STATIC ${C_SOURCES} ${C_WRAPPER_SOURCES} ${C_GENERATED_SOURCES} ${F_SOURCES})
target_compile_definitions(grackle PUBLIC LINUX H5_USE_16_API OMIT_LEGACY_INTERNAL_GRACKLE_FUNC)
target_include_directories(grackle PRIVATE ${C_INCLUDES} ${C_INCLUDES_DIRECTORY})
