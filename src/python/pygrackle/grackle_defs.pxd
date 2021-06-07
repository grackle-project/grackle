cdef extern from "grackle_types.h":
    # This does not need to be exactly correct, only of the right basic type
    ctypedef float gr_float

cdef extern from "grackle_chemistry_data.h":
    ctypedef struct c_chemistry_data "chemistry_data":
        int use_grackle
        int with_radiative_cooling
        int primordial_chemistry
        int dust_chemistry
        int metal_cooling
        int UVbackground
        char *grackle_data_file
        int cmb_temperature_floor
        double Gamma
        int h2_on_dust
        int use_dust_density_field
        int photoelectric_heating
        double photoelectric_heating_rate
        int use_isrf_field
        double interstellar_radiation_field
        int use_volumetric_heating_rate
        int use_specific_heating_rate
        int three_body_rate
        int cie_cooling
        int h2_optical_depth_approximation
        int ih2co
        int ipiht
        int h2_charge_exchange_rate
        int h2dust_rate
        int h2_h_cooling_rate
        int collisional_excitation_rates
        int collisional_ionisation_rates
        int recombination_cooling_rates
        int bremsstrahlung_cooling_rates
        double HydrogenFractionByMass
        double DeuteriumToHydrogenRatio
        double SolarMetalFractionByMass
        double local_dust_to_gas_ratio
        int NumberOfTemperatureBins
        int CaseBRecombination
        double TemperatureStart
        double TemperatureEnd
        int NumberOfDustTemperatureBins
        double DustTemperatureStart
        double DustTemperatureEnd
        int Compton_xray_heating
        int LWbackground_sawtooth_suppression
        double LWbackground_intensity
        double UVbackground_redshift_on
        double UVbackground_redshift_off
        double UVbackground_redshift_fullon
        double UVbackground_redshift_drop
        double cloudy_electron_fraction_factor
        int use_radiative_transfer
        int radiative_transfer_coupled_rate_solver
        int radiative_transfer_intermediate_step
        int radiative_transfer_hydrogen_only
        int self_shielding_method
        int H2_self_shielding

    ctypedef struct c_chemistry_data_storage "chemistry_data_storage":
        double *k1
        double *k2
        double *k3
        double *k4
        double *k5
        double *k6

        double *k7
        double *k8
        double *k9
        double *k10
        double *k11
        double *k12
        double *k13
        double *k14
        double *k15
        double *k16
        double *k17
        double *k18
        double *k19
        double *k20
        double *k21
        double *k22
        double *k23
        double *k13dd

        double k24
        double k25
        double k26

        double k27
        double k28
        double k29
        double k30
        double k31

        double *k50
        double *k51
        double *k52
        double *k53
        double *k54
        double *k55
        double *k56

        double *k57
        double *k58

        double *h2dust

        double *n_cr_n
        double *n_cr_d1
        double *n_cr_d2

        double *ceHI
        double *ceHeI
        double *ceHeII
        double *ciHI
        double *ciHeI
        double *ciHeIS
        double *ciHeII
        double *reHII
        double *reHeII1
        double *reHeII2
        double *reHeIII
        double *brem
        double comp
        double comp_xray
        double temp_xray

        double piHI
        double piHeI
        double piHeII

        double crsHI
        double crsHeI
        double crsHeII

        double *hyd01k
        double *h2k01
        double *vibh
        double *roth
        double *rotl
        double *GP99LowDensityLimit
        double *GP99HighDensityLimit

        double *GAHI
        double *GAH2
        double *GAHe
        double *GAHp
        double *GAel

        double *H2LTE

        double *HDlte
        double *HDlow

        double *cieco

        double gammah

        double *regr

        double gamma_isrf

        double *gas_grain
                                
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

    ctypedef struct c_field_data "grackle_field_data":
      int grid_rank;
      int *grid_dimension;
      int *grid_start;
      int *grid_end;
      gr_float grid_dx;
      gr_float *density;
      gr_float *HI_density;
      gr_float *HII_density;
      gr_float *HM_density;
      gr_float *HeI_density;
      gr_float *HeII_density;
      gr_float *HeIII_density;
      gr_float *H2I_density;
      gr_float *H2II_density;
      gr_float *DI_density;
      gr_float *DII_density;
      gr_float *HDI_density;
      gr_float *e_density;
      gr_float *metal_density;
      gr_float *dust_density;
      gr_float *internal_energy;
      gr_float *x_velocity;
      gr_float *y_velocity;
      gr_float *z_velocity;
      gr_float *volumetric_heating_rate;
      gr_float *specific_heating_rate;
      gr_float *RT_heating_rate;
      gr_float *RT_HI_ionization_rate;
      gr_float *RT_HeI_ionization_rate;
      gr_float *RT_HeII_ionization_rate;
      gr_float *RT_H2_dissociation_rate;
      gr_float *H2_self_shielding_length;
      gr_float *isrf_habing;

cdef extern from "grackle.h":
    c_chemistry_data _set_default_chemistry_parameters()

    int _initialize_chemistry_data(c_chemistry_data *my_chemistry,
                                   c_chemistry_data_storage *my_rates,
                                   c_code_units *my_units)

    int c_local_solve_chemistry "local_solve_chemistry"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                double dt_value)

    int c_local_calculate_cooling_time "local_calculate_cooling_time"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                gr_float *cooling_time)

    int c_local_calculate_gamma "local_calculate_gamma"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                gr_float *gamma)

    int c_local_calculate_pressure "local_calculate_pressure"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                gr_float *pressure)

    int c_local_calculate_temperature "local_calculate_temperature"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                gr_float *temperature)

    int c_local_calculate_dust_temperature "local_calculate_dust_temperature"(
                c_chemistry_data *my_chemistry,
                c_chemistry_data_storage *my_rates,
                c_code_units *my_units,
                c_field_data *my_fields,
                gr_float *dust_temperature)
