cdef extern from "grackle_types.h":
    # This does not need to be exactly correct, only of the right basic type
    ctypedef float gr_float

cdef extern from "grackle_chemistry_data.h":
    ctypedef struct c_chemistry_data "chemistry_data":
        double Gamma
        int use_grackle
        int with_radiative_cooling
        int primordial_chemistry
        int metal_cooling
        int h2_on_dust
        int cmb_temperature_floor
        char *grackle_data_file
        int three_body_rate
        int cie_cooling
        int h2_optical_depth_approximation
        int photoelectric_heating
        int CaseBRecombination
        int UVbackground
        double SolarMetalFractionByMass
        int use_volumetric_heating_rate
        int use_specific_heating_rate
        int use_radiative_transfer
        int self_shielding_method
        int H2_self_shielding

    ctypedef struct c_chemistry_data_storage "chemistry_data_storage":
        double k24
        double k25
        double k26
        double k27
        double k28
        double k29
        double k30
        double k31
        double hi_avg_crs
        double hei_avg_crs
        double heii_avg_crs

cdef extern from "grackle_types.h":
    ctypedef struct c_code_units "code_units":
      int comoving_coordinates
      double density_units
      double length_units
      double velocity_units
      double time_units
      double a_units
      double a_value

cdef extern from "grackle.h":
    c_chemistry_data _set_default_chemistry_parameters()

    int _initialize_chemistry_data(c_chemistry_data *my_chemistry,
                                   c_chemistry_data_storage *my_rates,
                                   c_code_units *my_units)

    int c_solve_chemistry "_solve_chemistry"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                double dt_value,
                double grid_dx,
                int grid_rank,
                int *grid_dimension,
                int *grid_start,
                int *grid_end,
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
                gr_float *volumetric_heating_rate,
                gr_float *specific_heating_rate,
                gr_float *RT_heating_rate,
                gr_float *RT_HI_ionization_rate,
                gr_float *RT_HeI_ionization_rate,
                gr_float *RT_HeII_ionization_rate,
                gr_float *RT_H2_dissociation_rate)

    int c_calculate_cooling_time "_calculate_cooling_time"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                int grid_rank,
                int *grid_dimension,
                int *grid_start,
                int *grid_end,
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
                gr_float *cooling_time,
                gr_float *RT_heating_rate,
                gr_float *volumetric_heating_rate,
                gr_float *specific_heating_rate)

    int c_calculate_gamma "_calculate_gamma"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                int grid_rank,
                int *grid_dimension,
                int *grid_start,
                int *grid_end,
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

    int c_calculate_pressure "_calculate_pressure"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                int grid_rank,
                int *grid_dimension,
                int *grid_start,
                int *grid_end,
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

    int c_calculate_temperature "_calculate_temperature"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                int grid_rank,
                int *grid_dimension,
                int *grid_start,
                int *grid_end,
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
