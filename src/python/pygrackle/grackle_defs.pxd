cdef extern from "grackle_macros.h":
    # This does not need to be exactly correct, only of the right basic type
    ctypedef float gr_float
    ctypedef int gr_int

cdef extern from "chemistry_data.h":
    ctypedef struct c_chemistry_data "chemistry_data":
        gr_float Gamma
        gr_int use_chemistry
        gr_int with_radiative_cooling
        gr_int primordial_chemistry
        gr_int metal_cooling
        gr_int h2_on_dust
        gr_int cmb_temperature_floor
        gr_int include_metal_heating
        char *grackle_data_file
        gr_int three_body_rate
        gr_int cie_cooling
        gr_int h2_optical_depth_approximation
        gr_int photoelectric_heating
        gr_int UVbackground
        # Most of the rest are not user-settable

cdef extern from "code_units.h":
    ctypedef struct c_code_units "code_units":
      gr_int comoving_coordinates
      gr_float density_units
      gr_float length_units
      gr_float time_units
      gr_float a_units

cdef extern from "grackle.h":
    c_chemistry_data set_default_chemistry_parameters()

    gr_int initialize_chemistry_data(c_chemistry_data &my_chemistry,
                                  c_code_units &my_units, gr_float a_value)

    gr_int update_UVbackground_rates(c_chemistry_data &my_chemistry,
                                  c_code_units &my_units, gr_float a_value)

    gr_int c_solve_chemistry "solve_chemistry"(
                c_chemistry_data &my_chemistry,
                c_code_units &my_units,
                gr_float a_value,
                gr_float dt_value,
                gr_int grid_rank,
                gr_int *grid_dimension,
                gr_int *grid_start,
                gr_int *grid_end,
                gr_float *density,
                gr_float *internal_energy,
                gr_float *x_velocity,
                gr_float *y_velocity,
                gr_float *z_velocity,
                gr_float *HI_density,
                gr_float *HII_density,
                gr_float *HM_density,
                gr_float *HeI_density,
                gr_float *HeII_density,
                gr_float *HeIII_density,
                gr_float *H2I_density,
                gr_float *H2II_density,
                gr_float *DI_density,
                gr_float *DII_density,
                gr_float *HDI_density,
                gr_float *e_density,
                gr_float *metal_density)

    gr_int c_calculate_cooling_time "calculate_cooling_time"(
                c_chemistry_data &my_chemistry,
                c_code_units &my_units,
                gr_float a_value,
                gr_int grid_rank,
                gr_int *grid_dimension,
                gr_int *grid_start,
                gr_int *grid_end,
                gr_float *density,
                gr_float *internal_energy,
                gr_float *x_velocity,
                gr_float *y_velocity,
                gr_float *z_velocity,
                gr_float *HI_density,
                gr_float *HII_density,
                gr_float *HM_density,
                gr_float *HeI_density,
                gr_float *HeII_density,
                gr_float *HeIII_density,
                gr_float *H2I_density,
                gr_float *H2II_density,
                gr_float *DI_density,
                gr_float *DII_density,
                gr_float *HDI_density,
                gr_float *e_density,
                gr_float *metal_density,
                gr_float *cooling_time)

    gr_int c_calculate_gamma "calculate_gamma"(
                c_chemistry_data &my_chemistry,
                c_code_units &my_units,
                gr_int grid_rank,
                gr_int *grid_dimension,
                gr_float *density,
                gr_float *internal_energy,
                gr_float *HI_density,
                gr_float *HII_density,
                gr_float *HM_density,
                gr_float *HeI_density,
                gr_float *HeII_density,
                gr_float *HeIII_density,
                gr_float *H2I_density,
                gr_float *H2II_density,
                gr_float *DI_density,
                gr_float *DII_density,
                gr_float *HDI_density,
                gr_float *e_density,
                gr_float *metal_density,
                gr_float *my_gamma)

    gr_int c_calculate_pressure "calculate_pressure"(
                c_chemistry_data &my_chemistry,
                c_code_units &my_units,
                gr_int grid_rank,
                gr_int *grid_dimension,
                gr_float *density,
                gr_float *internal_energy,
                gr_float *HI_density,
                gr_float *HII_density,
                gr_float *HM_density,
                gr_float *HeI_density,
                gr_float *HeII_density,
                gr_float *HeIII_density,
                gr_float *H2I_density,
                gr_float *H2II_density,
                gr_float *DI_density,
                gr_float *DII_density,
                gr_float *HDI_density,
                gr_float *e_density,
                gr_float *metal_density,
                gr_float *pressure)

    gr_int c_calculate_temperature "calculate_temperature"(
                c_chemistry_data &my_chemistry,
                c_code_units &my_units,
                gr_int grid_rank,
                gr_int *grid_dimension,
                gr_float *density,
                gr_float *internal_energy,
                gr_float *HI_density,
                gr_float *HII_density,
                gr_float *HM_density,
                gr_float *HeI_density,
                gr_float *HeII_density,
                gr_float *HeIII_density,
                gr_float *H2I_density,
                gr_float *H2II_density,
                gr_float *DI_density,
                gr_float *DII_density,
                gr_float *HDI_density,
                gr_float *e_density,
                gr_float *metal_density,
                gr_float *temperature)
