// See LICENSE file for license and copyright information

/// @file calc_tdust_3d.C
/// @brief Declares signature of calc_tdust_3d_g

// This file was initially generated automatically during conversion of the
// calc_tdust_3d_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "calc_tdust_3d.h"
#include "dust_props.hpp"
#include "fortran_func_wrappers.hpp"
#include "grackle.h"
#include "index_helper.h"
#include "internal_types.hpp"
#include "scale_fields.hpp"
#include "utils-cpp.hpp"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void calc_tdust_3d_g(
  gr_float* gas_temp_data_, gr_float* dust_temp_data_, int imetal,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, InternalGrUnits internalu
)
{
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  grackle::impl::View<gr_float***> metal(my_fields->metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(my_fields->dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_loc(my_fields->local_ISM_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C13(my_fields->ccsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C20(my_fields->ccsn20_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C25(my_fields->ccsn25_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_C30(my_fields->ccsn30_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F13(my_fields->fsn13_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F15(my_fields->fsn15_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F50(my_fields->fsn50_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_F80(my_fields->fsn80_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P170(my_fields->pisn170_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_P200(my_fields->pisn200_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_Y19(my_fields->y19_metal_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiM(my_fields->SiM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeM(my_fields->FeM_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg2SiO4(my_fields->Mg2SiO4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgSiO3(my_fields->MgSiO3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe3O4(my_fields->Fe3O4_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> AC(my_fields->AC_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2D(my_fields->SiO2_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> MgO(my_fields->MgO_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> FeS(my_fields->FeS_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al2O3(my_fields->Al2O3_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> reforg(my_fields->ref_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> volorg(my_fields->vol_org_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2Oice(my_fields->H2O_ice_dust_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const double mh_local_var = mh_grflt;

  // Set log values of start and end of lookup tables

  const double logtem0  = std::log(my_chemistry->TemperatureStart);
  const double logtem9  = std::log(my_chemistry->TemperatureEnd);
  const double dlogtem  = (logtem9 - logtem0)/(double)(my_chemistry->NumberOfTemperatureBins-1);

  // Set unit-related quantities
  const double dom = internalu_calc_dom_(internalu);

  // Set CMB temperature
  const double zr = (gr_float)(1.)/(internalu.a_value*internalu.a_units) - (gr_float)(1.);
  const double trad = (gr_float)(2.73) * ((gr_float)(1.) + zr);

  // Loop over zones, and do an entire i-column in one go
  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  // Convert densities to 'proper' from comoving
  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(1.0)/(gr_float)std::pow(internalu.a_value,3);
    grackle::impl::scale_fields_dust(my_chemistry, my_fields, imetal, factor);
  }

  OMP_PRAGMA("omp parallel")
  {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

    grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> gas_temp(gas_temp_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> dust_temp(dust_temp_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> isrf_habing(my_fields->isrf_habing, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    grackle::impl::View<gr_float***> SiM_temp(my_fields->SiM_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> FeM_temp(my_fields->FeM_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> Mg2SiO4_temp(my_fields->Mg2SiO4_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> MgSiO3_temp(my_fields->MgSiO3_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> Fe3O4_temp(my_fields->Fe3O4_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> AC_temp(my_fields->AC_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> SiO2D_temp(my_fields->SiO2_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> MgO_temp(my_fields->MgO_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> FeS_temp(my_fields->FeS_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> Al2O3_temp(my_fields->Al2O3_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> reforg_temp(my_fields->ref_org_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> volorg_temp(my_fields->vol_org_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> H2Oice_temp(my_fields->H2O_ice_dust_temperature, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    // buffers that store values computed within an index-range
    std::vector<double> metallicity(my_fields->grid_dimension[0]);
    std::vector<double> dust2gas(my_fields->grid_dimension[0]);
    std::vector<double> tgas(my_fields->grid_dimension[0]);
    std::vector<double> tdust(my_fields->grid_dimension[0]);
    std::vector<double> nh(my_fields->grid_dimension[0]);
    std::vector<double> gasgr(my_fields->grid_dimension[0]);
    std::vector<double> gasgr_tdust(my_fields->grid_dimension[0]);
    std::vector<double> myisrf(my_fields->grid_dimension[0]);
    std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    grackle::impl::InternalDustPropBuf internal_dust_prop_buf =
      grackle::impl::new_InternalDustPropBuf(
        my_fields->grid_dimension[0], my_rates->gr_N[1]
      );

    // these next buffer variables hold values that are computed as side-effect
    // (we don't really care about them in the current context)

    // opacity coefficients for each dust grain (the product of opacity
    // coefficient & gas mass density is the linear absortpion coefficient)
    grackle::impl::GrainSpeciesCollection grain_kappa =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    // closely related to grain_kappa
    std::vector<double> kappa_tot(my_fields->grid_dimension[0]);

    // holds the gas/grain-species heat transfer rates
    grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);


    // The following for-loop is a flattened loop over every k,j combination.
    // OpenMP divides this loop between all threads. Within the loop, we
    // complete calculations for the constructed index-range construct
    // (an index range corresponds to an "i-slice")
    OMP_PRAGMA("omp for schedule(runtime)")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);
      const int k = idx_range.k;
      const int j = idx_range.j;


      // Set itmask to true for entire idx_range
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask_metal[i] = MASK_TRUE;
      }

      // Set itmask to false for metal-poor cells
      if (imetal == 1) {
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (metal(i,j,k) < 1.e-9 * d(i,j,k))  {
            itmask_metal[i] = MASK_FALSE;
          }
        }
      }

      // Compute grain size increment

      if ( (my_chemistry->use_dust_density_field > 0)  &&  (my_chemistry->dust_species > 0) )  {

        f_wrap::calc_grain_size_increment_1d (
          dom, idx_range, itmask_metal.data(), my_chemistry, my_rates,
          my_fields, internal_dust_prop_buf
        );

      }

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if(itmask_metal[i] != MASK_FALSE)  {
          // Calculate metallicity

          if (imetal == 1)  {
            metallicity[i] = metal(i,j,k) / d(i,j,k) / my_chemistry->SolarMetalFractionByMass;
          }

          // Calculate dust to gas ratio

          //       if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
          //        if (idspecies .gt. 0) then
          //           dust(i,j,k) = MgSiO3  (i,j,k)
          // &                     + AC      (i,j,k)
          //        endif
          //        if (idspecies .gt. 1) then
          //           dust(i,j,k) = dust(i,j,k)
          // &                     + SiM     (i,j,k)
          // &                     + FeM     (i,j,k)
          // &                     + Mg2SiO4 (i,j,k)
          // &                     + Fe3O4   (i,j,k)
          // &                     + SiO2D   (i,j,k)
          // &                     + MgO     (i,j,k)
          // &                     + FeS     (i,j,k)
          // &                     + Al2O3   (i,j,k)
          //        endif
          //        if (idspecies .gt. 2) then
          //           dust(i,j,k) = dust(i,j,k)
          // &                     + reforg  (i,j,k)
          // &                     + volorg  (i,j,k)
          // &                     + H2Oice  (i,j,k)
          //        endif
          //       endif

          if (my_chemistry->use_dust_density_field > 0)  {
            dust2gas[i] = dust(i,j,k) / d(i,j,k);
          } else {
            dust2gas[i] = my_chemistry->local_dust_to_gas_ratio * metallicity[i];
          }

          // Compute interstellar radiation field

          if (my_chemistry->use_isrf_field > 0)  {
            myisrf[i] = isrf_habing(i,j,k);
          } else {
            myisrf[i] = my_chemistry->interstellar_radiation_field;
          }


          // Compute hydrogen number density

          nh[i] = HI(i,j,k) + HII(i,j,k);
          if (my_chemistry->primordial_chemistry > 1)  {
            nh[i] = nh[i] + H2I(i,j,k) + H2II(i,j,k);
          }

          // We have not converted to proper, so use urho and not dom

          nh[i] = nh[i] * internalu.urho / mh_local_var;

          // Compute log temperature and truncate if above/below table max/min

          tgas[i]   = gas_temp(i,j,k);
          logTlininterp_buf.logtem[i] = std::log(tgas[i]);
          logTlininterp_buf.logtem[i] = std::fmax(logTlininterp_buf.logtem[i], logtem0);
          logTlininterp_buf.logtem[i] = std::fmin(logTlininterp_buf.logtem[i], logtem9);

          // Compute index into the table and precompute parts of linear interp

          logTlininterp_buf.indixe[i] = std::fmin(
            my_chemistry->NumberOfTemperatureBins-1,
            std::fmax(1,(long long)((logTlininterp_buf.logtem[i]-logtem0)/dlogtem)+1)
          );
          logTlininterp_buf.t1[i] = (logtem0 + (logTlininterp_buf.indixe[i] - 1)*dlogtem);
          logTlininterp_buf.t2[i] = (logtem0 + (logTlininterp_buf.indixe[i]    )*dlogtem);
          logTlininterp_buf.tdef[i] = (logTlininterp_buf.logtem[i] - logTlininterp_buf.t1[i]) / (logTlininterp_buf.t2[i] - logTlininterp_buf.t1[i]);

        }
      }

      // Compute dust temperature(s) in the index-range
      f_wrap::calc_all_tdust_gasgr_1d_g(
        trad, tgas.data(), tdust.data(), metallicity.data(), dust2gas.data(),
        nh.data(), gasgr_tdust.data(), itmask_metal.data(), internalu.coolunit,
        gasgr.data(), myisrf.data(), kappa_tot.data(), my_chemistry, my_rates,
        my_fields, idx_range, grain_temperatures, gas_grainsp_heatrate,
        grain_kappa, logTlininterp_buf, internal_dust_prop_buf
      );

      // Copy slice values back to grid

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        if (itmask_metal[i] != MASK_FALSE) {
          if (my_chemistry->use_multiple_dust_temperatures == 0)  {
            dust_temp(i,j,k) = tdust[i];
          } else {
            if (my_chemistry->dust_species > 0)  {
              MgSiO3_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i];
              AC_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i];
            }
            if (my_chemistry->dust_species > 1)  {
              SiM_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i];
              FeM_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i];
              Mg2SiO4_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i];
              Fe3O4_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i];
              SiO2D_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i];
              MgO_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i];
              FeS_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i];
              Al2O3_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i];
            }
            if (my_chemistry->dust_species > 2)  {
              reforg_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][i];
              volorg_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i];
              H2Oice_temp(i,j,k) = grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i];
            }
          }
        }
      }

    }

    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
    grackle::impl::drop_InternalDustPropBuf(&internal_dust_prop_buf);
    grackle::impl::drop_GrainSpeciesCollection(&grain_kappa);
    grackle::impl::drop_GrainSpeciesCollection(&gas_grainsp_heatrate);
  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from 'proper'
  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)std::pow(internalu.a_value,3);
    grackle::impl::scale_fields_dust(my_chemistry, my_fields, imetal, factor);
  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
