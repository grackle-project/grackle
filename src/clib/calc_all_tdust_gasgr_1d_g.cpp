//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_all_tdust_gasgr_1d_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_all_tdust_gasgr_1d_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "calc_all_tdust_gasgr_1d_g.hpp"

void grackle::impl::calc_all_tdust_gasgr_1d_g(
  double trad, double* tgas, double* tdust, const double* metallicity,
  const double* dust2gas, double* nh, double* gasgr_tdust, gr_mask_type* itmask_metal,
  double coolunit, double* gasgr, double* myisrf, double* kptot,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf,
  grackle::impl::GrainSpeciesCollection grain_kappa
)
{
  const double mh_local_var = mh_grflt;
  int i;
  
  // Cooling/heating slice locals

  std::vector<double> gasgr_tSiM(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tFeM(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tAC(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tMgO(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tFeS(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_treforg(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tvolorg(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tH2Oice(my_fields->grid_dimension[0]);
  double fv2k, fac;
  std::vector<double> mygisrf(my_fields->grid_dimension[0]);
  std::vector<double> gisrfSiM(my_fields->grid_dimension[0]);
  std::vector<double> gisrfFeM(my_fields->grid_dimension[0]);
  std::vector<double> gisrfMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> gisrfMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> gisrfFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> gisrfAC(my_fields->grid_dimension[0]);
  std::vector<double> gisrfSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> gisrfMgO(my_fields->grid_dimension[0]);
  std::vector<double> gisrfFeS(my_fields->grid_dimension[0]);
  std::vector<double> gisrfAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> gisrfreforg(my_fields->grid_dimension[0]);
  std::vector<double> gisrfvolorg(my_fields->grid_dimension[0]);
  std::vector<double> gisrfH2Oice(my_fields->grid_dimension[0]);

  // Calculate heating from interstellar radiation field
  //  -> this is ONLY used when `itmask_metal .eq. MASK_TRUE`

  for (i = idx_range.i_start; i<=idx_range.i_end; i++) {
    if ( itmask_metal[i] != MASK_FALSE )  {

      if (my_chemistry->dust_species == 0 )  {
        if (my_chemistry->use_dust_density_field > 0)  {
          mygisrf[i] = my_rates->gamma_isrf
                     * my_chemistry->local_dust_to_gas_ratio / dust2gas[i] * metallicity[i];
          // ! correct with the depletion or enhancement of condensation rate.
        } else {
          mygisrf[i] = my_rates->gamma_isrf;
        }

      } else {

        if (my_chemistry->use_multiple_dust_temperatures == 0)  {
            
          mygisrf[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.sigma_per_gas_mass_tot[i];
          // ! in UV, absorption coefficient Q ~ 1 (Goldsmith 2001)
          // ! so we use the geometrical cross-section of grains [cgs]

        } else {

          if (my_chemistry->dust_species > 0)  {
            gisrfMgSiO3[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust][i];
            gisrfAC[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust][i];
          }
          if (my_chemistry->dust_species > 1)  {
            gisrfSiM[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust][i];
            gisrfFeM[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust][i];
            gisrfMg2SiO4[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust][i];
            gisrfFe3O4[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust][i];
            gisrfSiO2D[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust][i];
            gisrfMgO[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust][i];
            gisrfFeS[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust][i];
            gisrfAl2O3[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust][i];
          }
          if (my_chemistry->dust_species > 2)  {
            gisrfreforg[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust][i];
            gisrfvolorg[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust][i];
            gisrfH2Oice[i] = my_rates->gamma_isrf2 * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust][i];
          }

        }

      }

    }
  }

  // --- Gas to grain heat transfer ---

  // Look up gas/grain heat transfer rates

  for (i = idx_range.i_start; i<=idx_range.i_end; i++) {
    if ( itmask_metal[i] != MASK_FALSE )  {

      if(my_chemistry->dust_species == 0)  {

        gasgr[i] = my_rates->gas_grain[logTlininterp_buf.indixe[i]-1] + logTlininterp_buf.tdef[i]
             *(my_rates->gas_grain[logTlininterp_buf.indixe[i]+1-1] -my_rates->gas_grain[logTlininterp_buf.indixe[i]-1]);

        gasgr_tdust[i] = (dust2gas[i] / metallicity[i])
                             * gasgr[i] * coolunit / mh_local_var;

      } else {

        fv2k = my_rates->gas_grain2[logTlininterp_buf.indixe[i]-1] + logTlininterp_buf.tdef[i]
             *(my_rates->gas_grain2[logTlininterp_buf.indixe[i]+1-1] -my_rates->gas_grain2[logTlininterp_buf.indixe[i]-1]);

        fac = coolunit / mh_local_var;

        if ( my_chemistry->use_multiple_dust_temperatures == 0 )  {

          gasgr[i] = fv2k * internal_dust_prop_buf.sigma_per_gas_mass_tot[i];

          gasgr_tdust[i] = gasgr[i] * fac;
                 
        } else {

          if (my_chemistry->dust_species > 0)  {
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust][i];
          }
          if (my_chemistry->dust_species > 1)  {
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust][i];
          }
          if (my_chemistry->dust_species > 2)  {
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust][i];
            gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] = fv2k * internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust][i];
          }

          if (my_chemistry->dust_species > 0)  {
            gasgr_tMgSiO3[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust][i] * fac;
            gasgr_tAC[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust][i] * fac;
          }
          if (my_chemistry->dust_species > 1)  {
            gasgr_tSiM[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust][i] * fac;
            gasgr_tFeM[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust][i] * fac;
            gasgr_tMg2SiO4[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust][i] * fac;
            gasgr_tFe3O4[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust][i] * fac;
            gasgr_tSiO2D[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust][i] * fac;
            gasgr_tMgO[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust][i] * fac;
            gasgr_tFeS[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust][i] * fac;
            gasgr_tAl2O3[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust][i] * fac;
          }
          if (my_chemistry->dust_species > 2)  {
            gasgr_treforg[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust][i] * fac;
            gasgr_tvolorg[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust][i] * fac;
            gasgr_tH2Oice[i] = gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust][i] * fac;
          }

        }

      }

    }
  }

  // Compute dust temperature

  if (my_chemistry->use_multiple_dust_temperatures == 0)  {

     FORTRAN_NAME(calc_tdust_1d_g)(tdust, tgas, nh, gasgr_tdust,
      mygisrf.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
      my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.dyntab_kappa_tot, kptot, &my_chemistry->dust_species);

  } else {

    if (my_chemistry->dust_species > 0)  {
      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], tgas, nh, gasgr_tMgSiO3.data(),
        gisrfMgSiO3.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgSiO3_dust],
        grain_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::AC_dust], tgas, nh, gasgr_tAC.data(),
        gisrfAC.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust],
        grain_kappa.data[OnlyGrainSpLUT::AC_dust], &my_chemistry->dust_species);
    }

    if (my_chemistry->dust_species > 1)  {
      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], tgas, nh, gasgr_tSiM.data(),
        gisrfSiM.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust],
        grain_kappa.data[OnlyGrainSpLUT::SiM_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], tgas, nh, gasgr_tFeM.data(),
        gisrfFeM.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust],
        grain_kappa.data[OnlyGrainSpLUT::FeM_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], tgas, nh, gasgr_tMg2SiO4.data(),
        gisrfMg2SiO4.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust],
        grain_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust], tgas, nh, gasgr_tFe3O4.data(),
        gisrfFe3O4.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Fe3O4_dust],
        grain_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], tgas, nh, gasgr_tSiO2D.data(),
        gisrfSiO2D.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust],
        grain_kappa.data[OnlyGrainSpLUT::SiO2_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], tgas, nh, gasgr_tMgO.data(),
        gisrfMgO.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust],
        grain_kappa.data[OnlyGrainSpLUT::MgO_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], tgas, nh, gasgr_tFeS.data(),
        gisrfFeS.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust],
        grain_kappa.data[OnlyGrainSpLUT::FeS_dust], &my_chemistry->dust_species);

      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], tgas, nh, gasgr_tAl2O3.data(),
        gisrfAl2O3.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Al2O3_dust],
        grain_kappa.data[OnlyGrainSpLUT::Al2O3_dust], &my_chemistry->dust_species);
    }

    if (my_chemistry->dust_species > 2)  {
      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust], tgas, nh, gasgr_treforg.data(),
        gisrfreforg.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::ref_org_dust],
        grain_kappa.data[OnlyGrainSpLUT::ref_org_dust], &my_chemistry->dust_species);
          
      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], tgas, nh, gasgr_tvolorg.data(),
        gisrfvolorg.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::vol_org_dust],
        grain_kappa.data[OnlyGrainSpLUT::vol_org_dust], &my_chemistry->dust_species);
          
      FORTRAN_NAME(calc_tdust_1d_g)(grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], tgas, nh, gasgr_tH2Oice.data(),
        gisrfH2Oice.data(), myisrf, itmask_metal, &trad, &my_fields->grid_dimension[0], &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1,
        my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td, internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::H2O_ice_dust],
        grain_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], &my_chemistry->dust_species);
    }

  }
}