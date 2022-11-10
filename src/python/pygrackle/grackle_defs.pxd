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

        int metal_chemistry

        int grain_growth

        int multi_metals

        int metal_abundances

        int dust_species

        int dust_temperature_multi

        int dust_sublimation

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
        double HydrogenFractionByMass
        double DeuteriumToHydrogenRatio
        double SolarMetalFractionByMass
        double local_dust_to_gas_ratio

        int SN0_N
        double *SN0_XC 
        double *SN0_XO 
        double *SN0_XMg
        double *SN0_XAl
        double *SN0_XSi
        double *SN0_XS 
        double *SN0_XFe
        double *SN0_fC 
        double *SN0_fO 
        double *SN0_fMg
        double *SN0_fAl
        double *SN0_fSi
        double *SN0_fS 
        double *SN0_fFe
        double *SN0_fSiM
        double *SN0_fFeM
        double *SN0_fMg2SiO4
        double *SN0_fMgSiO3
        double *SN0_fFe3O4
        double *SN0_fAC
        double*SN0_fSiO2D
        double*SN0_fMgO
        double*SN0_fFeS
        double*SN0_fAl2O3
        double *SN0_freforg 
        double*SN0_fvolorg 
        double*SN0_fH2Oice
        double *SN0_r0SiM
        double *SN0_r0FeM
        double *SN0_r0Mg2SiO4
        double *SN0_r0MgSiO3
        double *SN0_r0Fe3O4
        double *SN0_r0AC
        double*SN0_r0SiO2D
        double*SN0_r0MgO
        double*SN0_r0FeS
        double*SN0_r0Al2O3
        double *SN0_r0reforg 
        double*SN0_r0volorg 
        double*SN0_r0H2Oice

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
        int radiative_transfer_H2II_diss
        int radiative_transfer_HDI_diss
        int radiative_transfer_metal_ion
        int radiative_transfer_metal_diss

        int radiative_transfer_use_H2_shielding

        int self_shielding_method

        int H2_self_shielding

        int h2_charge_exchange_rate

        int h2_dust_rate

        int h2_h_cooling_rate

        int collisional_excitation_rates
        int collisional_ionisation_rates
        int recombination_cooling_rates
        int bremsstrahlung_cooling_rates

        int use_palla_salpeter_stahler_1983
        int use_stancil_lepp_dalgarno_1998
        int use_omukai_gas_grain
        int use_uniform_grain_dist_gamma_isrf

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

        double *k125
        double *k129
        double *k130
        double *k131
        double *k132
        double *k133
        double *k134
        double *k135
        double *k136
        double *k137
        double *k148
        double *k149
        double *k150
        double *k151
        double *k152
        double *k153

        double *kz15
        double *kz16
        double *kz17
        double *kz18
        double *kz19
        double *kz20
        double *kz21
        double *kz22
        double *kz23
        double *kz24
        double *kz25
        double *kz26
        double *kz27
        double *kz28
        double *kz29
        double *kz30
        double *kz31
        double *kz32
        double *kz33
        double *kz34
        double *kz35
        double *kz36
        double *kz37
        double *kz38
        double *kz39
        double *kz40
        double *kz41
        double *kz42
        double *kz43
        double *kz44
        double *kz45
        double *kz46
        double *kz47
        double *kz48
        double *kz49
        double *kz50
        double *kz51
        double *kz52
        double *kz53
        double *kz54

        double *h2dust
        double *h2dustS
        double *h2dustC

        double *grain_growth

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
        double gamma_isrf2

        double *gas_grain
        double *gas_grain2

        double *cieY06

        int *LH2_N
        int LH2_Size
        double *LH2_D
        double *LH2_T
        double *LH2_H
        double LH2_dD
        double LH2_dT
        double LH2_dH
        double *LH2_L
        int *LHD_N
        int LHD_Size
        double *LHD_D
        double *LHD_T
        double *LHD_H
        double LHD_dD
        double LHD_dT
        double LHD_dH
        double *LHD_L

        int *LCI_N
        int LCI_Size
        double *LCI_D
        double *LCI_T
        double *LCI_H
        double LCI_dD
        double LCI_dT
        double LCI_dH
        double *LCI_L
        int *LCII_N
        int LCII_Size
        double *LCII_D
        double *LCII_T
        double *LCII_H
        double LCII_dD
        double LCII_dT
        double LCII_dH
        double *LCII_L
        int *LOI_N
        int LOI_Size
        double *LOI_D
        double *LOI_T
        double *LOI_H
        double LOI_dD
        double LOI_dT
        double LOI_dH
        double *LOI_L

        int *LCO_N
        int LCO_Size
        double *LCO_D
        double *LCO_T
        double *LCO_H
        double LCO_dD
        double LCO_dT
        double LCO_dH
        double *LCO_L
        int *LOH_N
        int LOH_Size
        double *LOH_D
        double *LOH_T
        double *LOH_H
        double LOH_dD
        double LOH_dT
        double LOH_dH
        double *LOH_L
        int *LH2O_N
        int LH2O_Size
        double *LH2O_D
        double *LH2O_T
        double *LH2O_H
        double LH2O_dD
        double LH2O_dT
        double LH2O_dH
        double *LH2O_L

        int *alphap_N
        int alphap_Size
        double *alphap_D
        double *alphap_T
        double alphap_dD
        double alphap_dT
        double *alphap_Data

        int *grain_N
        int grain_Size
        double *grain_D
        double *grain_T
        double grain_dD
        double grain_dT
        double *Hgrai
        double *Tgrai
        double *Ograi
        double *Lgrain

        int *gr_N
        int gr_Size
        double gr_dT
        double *gr_Td
        double *SN0_kpSiM
        double *SN0_kpFeM
        double *SN0_kpMg2SiO4
        double *SN0_kpMgSiO3
        double *SN0_kpFe3O4
        double *SN0_kpAC
        double *SN0_kpSiO2D
        double *SN0_kpMgO
        double *SN0_kpFeS
        double *SN0_kpAl2O3
        double *SN0_kpreforg 
        double *SN0_kpvolorg 
        double *SN0_kpH2Oice
                                
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
