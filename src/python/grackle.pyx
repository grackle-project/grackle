cdef extern from "chemistry_data.h":
    ctypedef struct chemistry_data:
        float Gamma
        int use_chemistry
        int primordial_chemistry
        int metal_cooling
        int h2_on_dust
        int cmb_temperature_floor
        int include_metal_heating
        char *cloudy_table_file
        int three_body_rate
        int cie_cooling
        int h2_optical_depth_approximation
        int photoelectric_heating
        # Most of the rest are not user-settable

cdef extern from "code_units.h":
    ctypedef struct code_units:
      int comoving_coordinates
      float density_units
      float length_units
      float time_units
      float a_units

cdef extern from "grackle.h":
    chemistry_data set_default_chemistry_parameters()

    int initialize_chemistry_data(chemistry_data &my_chemistry,
                                  code_units &my_units, float a_value)

    int solve_chemistry(chemistry_data &my_chemistry,
                code_units &my_units,
                float a_value, float dt_value,
                int grid_rank, int *grid_dimension,
                int *grid_start, int *grid_end,
                float *density, float *internal_energy,
                float *x_velocity, float *y_velocity, float  *z_velocity,
                float *HI_density, float *HII_density, float *HM_density,
                float *HeI_density, float *HeII_density, float *HeIII_density,
                float *H2I_density, float *H2II_density,
                float *DI_density, float *DII_density, float *HDI_density,
                float *e_density, float *metal_density)

    int calculate_cooling_time(chemistry_data &my_chemistry,
                   code_units &my_units,
                   float a_value, float dt_value,
                   int grid_rank, int *grid_dimension,
                   int *grid_start, int *grid_end,
                   float *density, float *internal_energy,
                   float *x_velocity, float *y_velocity, float  *z_velocity,
                   float *HI_density, float *HII_density, float *HM_density,
                   float *HeI_density, float *HeII_density, float *HeIII_density,
                   float *H2I_density, float *H2II_density,
                   float *DI_density, float *DII_density, float *HDI_density,
                   float *e_density, float *metal_density,
                   float *cooling_time)

    int calculate_gamma(chemistry_data &my_chemistry,
                        code_units &my_units,
                        int grid_rank, int *grid_dimension,
                        float *density, float *internal_energy,
                        float *HI_density, float *HII_density, float *HM_density,
                        float *HeI_density, float *HeII_density, float *HeIII_density,
                        float *H2I_density, float *H2II_density,
                        float *DI_density, float *DII_density, float *HDI_density,
                        float *e_density, float *metal_density,
                        float *my_gamma)

    int calculate_pressure(chemistry_data &my_chemistry,
                           code_units &my_units,
                           int grid_rank, int *grid_dimension,
                           float *density, float *internal_energy,
                           float *HI_density, float *HII_density, float *HM_density,
                           float *HeI_density, float *HeII_density, float *HeIII_density,
                           float *H2I_density, float *H2II_density,
                           float *DI_density, float *DII_density, float *HDI_density,
                           float *e_density, float *metal_density,
                           float *pressure)

    int calculate_temperature(chemistry_data &my_chemistry,
                              code_units &my_units,
                              int grid_rank, int *grid_dimension,
                              float *density, float *internal_energy,
                              float *HI_density, float *HII_density, float *HM_density,
                              float *HeI_density, float *HeII_density, float *HeIII_density,
                              float *H2I_density, float *H2II_density,
                              float *DI_density, float *DII_density, float *HDI_density,
                              float *e_density, float *metal_density,
                              float *temperature)
