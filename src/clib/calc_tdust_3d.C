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

  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> gas_temp(gas_temp_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_temp(dust_temp_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> isrf_habing(my_fields->isrf_habing, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  double dom;
  std::vector<double> metallicity(my_fields->grid_dimension[0]);
  std::vector<double> dust2gas(my_fields->grid_dimension[0]);
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

  // Parameters

  const double mh_local_var = mh_grflt;

  // Locals

  double trad, zr, logtem0, logtem9, dlogtem, coolunit, dbase1, tbase1, xbase1;

  // Slice locals
 
  std::vector<double> tgas(my_fields->grid_dimension[0]);
  std::vector<double> tdust(my_fields->grid_dimension[0]);
  std::vector<double> nh(my_fields->grid_dimension[0]);
  std::vector<double> gasgr(my_fields->grid_dimension[0]);
  std::vector<double> gasgr_tdust(my_fields->grid_dimension[0]);
  std::vector<double> myisrf(my_fields->grid_dimension[0]);

  // Iteration mask for multi_cool

  std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  // Set log values of start and end of lookup tables

  logtem0  = std::log(my_chemistry->TemperatureStart);
  logtem9  = std::log(my_chemistry->TemperatureEnd);
  dlogtem  = (logtem9 - logtem0)/(double)(my_chemistry->NumberOfTemperatureBins-1);

  // Set units

  dom      = internalu.urho*(std::pow(internalu.a_value,3))/mh_local_var;
  tbase1   = internalu.tbase1;
  xbase1   = internalu.uxyz/(internalu.a_value*internalu.a_units);    // uxyz is [x]*a      = [x]*[a]*a'        '
  dbase1   = internalu.urho*std::pow((internalu.a_value*internalu.a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '
  coolunit = (std::pow(internalu.a_units,5) * std::pow(xbase1,2) * std::pow(mh_local_var,2)) / (std::pow(tbase1,3) * dbase1);
  zr       = (gr_float)(1.)/(internalu.a_value*internalu.a_units) - (gr_float)(1.);

  // Set CMB temperature

  trad = (gr_float)(2.73) * ((gr_float)(1.) + zr);

  // Loop over zones, and do an entire i-column in one go
  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  if (internalu.extfields_in_comoving == 1)  {

    //_// PORT: #ifdef _OPENMP
    //_// PORT: !$omp parallel do schedule(runtime)
    //_// PORT: #endif
    //_// TODO_USE: OMP_PRAGMA("omp parallel")
    {
      //_// TODO: move relevant variable declarations to here to replace OMP private

      //_// TODO_USE: OMP_PRAGMA("omp for")
      for (int t = 0; t < idx_helper.outer_ind_size; t++) {
        // construct an index-range corresponding to "i-slice"
        const IndexRange idx_range = make_idx_range_(t, &idx_helper);
        const int k = idx_range.k;
        const int j = idx_range.j;

        if (imetal == 1)  {
          for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
            metal(ip1-1,j,k) = metal(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
            if (my_chemistry->multi_metals > 0)  {
              metal_loc(ip1-1,j,k) = metal_loc(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C13(ip1-1,j,k) = metal_C13(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C20(ip1-1,j,k) = metal_C20(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C25(ip1-1,j,k) = metal_C25(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C30(ip1-1,j,k) = metal_C30(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F13(ip1-1,j,k) = metal_F13(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F15(ip1-1,j,k) = metal_F15(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F50(ip1-1,j,k) = metal_F50(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F80(ip1-1,j,k) = metal_F80(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_P170(ip1-1,j,k)=metal_P170(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_P200(ip1-1,j,k)=metal_P200(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_Y19(ip1-1,j,k) = metal_Y19(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
            }
          }
        }
        if (my_chemistry->use_dust_density_field == 1)  {
          for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
            dust(ip1-1,j,k) = dust(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              // !            if (metal(ip1,j,k) .gt. 1.d-9 * d(ip1,j,k)) then
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(ip1-1,j,k)  = MgSiO3(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                AC(ip1-1,j,k)      = AC(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(ip1-1,j,k)     = SiM(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                FeM(ip1-1,j,k)     = FeM(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                Mg2SiO4(ip1-1,j,k) = Mg2SiO4(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                Fe3O4(ip1-1,j,k)   = Fe3O4(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                SiO2D(ip1-1,j,k)   = SiO2D(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                MgO(ip1-1,j,k)     = MgO(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                FeS(ip1-1,j,k)     = FeS(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                Al2O3(ip1-1,j,k)   = Al2O3(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(ip1-1,j,k)  = reforg(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                volorg(ip1-1,j,k)  = volorg(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
                H2Oice(ip1-1,j,k)  = H2Oice(ip1-1,j,k)/(gr_float)(std::pow(internalu.a_value,3) );
              }
              // !            endif
            }
          }
        }
      }
    }  // OMP_PRAGMA("omp parallel")
 
  }

  // Loop over slices (in the k-direction)

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: !$omp parallel do schedule(runtime) private(
  //_// PORT: !$omp&   tgas, tdust, nh, gasgr, myisrf,
  //_// PORT: !$omp&   itmask )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    //_// TODO: move relevant variable declarations to here to replace OMP private
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


    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);
      const int k = idx_range.k;
      const int j = idx_range.j;

      for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {

        // Set itmask to all true

        itmask[ip1-1] = MASK_TRUE;

      }

      // Iteration mask for metal-rich cells

      if (imetal == 1)  {
        for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
          if (metal(ip1-1,j,k) < 1.e-9 * d(ip1-1,j,k))  {
            itmask[ip1-1] = MASK_FALSE;
          }
        }
      }

      // Compute grain size increment

      if ( (my_chemistry->use_dust_density_field > 0)  &&  (my_chemistry->dust_species > 0) )  {

        f_wrap::calc_grain_size_increment_1d (
          dom, idx_range, itmask.data(), my_chemistry, my_rates, my_fields,
          internal_dust_prop_buf
        );

      }

      for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
        if(itmask[ip1-1] != MASK_FALSE)  {
          // Calculate metallicity

          if (imetal == 1)  {
            metallicity[ip1-1] = metal(ip1-1,j,k) / d(ip1-1,j,k) / my_chemistry->SolarMetalFractionByMass;
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
            dust2gas[ip1-1] = dust(ip1-1,j,k) / d(ip1-1,j,k);
          } else {
            dust2gas[ip1-1] = my_chemistry->local_dust_to_gas_ratio * metallicity[ip1-1];
          }

          // Compute interstellar radiation field

          if (my_chemistry->use_isrf_field > 0)  {
            myisrf[ip1-1] = isrf_habing(ip1-1,j,k);
          } else {
            myisrf[ip1-1] = my_chemistry->interstellar_radiation_field;
          }


          // Compute hydrogen number density

          nh[ip1-1] = HI(ip1-1,j,k) + HII(ip1-1,j,k);
          if (my_chemistry->primordial_chemistry > 1)  {
            nh[ip1-1] = nh[ip1-1] + H2I(ip1-1,j,k) + H2II(ip1-1,j,k);
          }

          // We have not converted to proper, so use urho and not dom

          nh[ip1-1] = nh[ip1-1] * internalu.urho / mh_local_var;

          // Compute log temperature and truncate if above/below table max/min

          tgas[ip1-1]   = gas_temp(ip1-1,j,k);
          logTlininterp_buf.logtem[ip1-1] = std::log(tgas[ip1-1]);
          logTlininterp_buf.logtem[ip1-1] = std::fmax(logTlininterp_buf.logtem[ip1-1], logtem0);
          logTlininterp_buf.logtem[ip1-1] = std::fmin(logTlininterp_buf.logtem[ip1-1], logtem9);

          // Compute index into the table and precompute parts of linear interp

          logTlininterp_buf.indixe[ip1-1] = std::fmin(
            my_chemistry->NumberOfTemperatureBins-1,
            std::fmax(1,(long long)((logTlininterp_buf.logtem[ip1-1]-logtem0)/dlogtem)+1)
          );
          logTlininterp_buf.t1[ip1-1] = (logtem0 + (logTlininterp_buf.indixe[ip1-1] - 1)*dlogtem);
          logTlininterp_buf.t2[ip1-1] = (logtem0 + (logTlininterp_buf.indixe[ip1-1]    )*dlogtem);
          logTlininterp_buf.tdef[ip1-1] = (logTlininterp_buf.logtem[ip1-1] - logTlininterp_buf.t1[ip1-1]) / (logTlininterp_buf.t2[ip1-1] - logTlininterp_buf.t1[ip1-1]);

        }
      }

      // Compute dust temperature(s) in the index-range
      f_wrap::calc_all_tdust_gasgr_1d_g(
        trad, tgas.data(), tdust.data(), metallicity.data(), dust2gas.data(),
        nh.data(), gasgr_tdust.data(), itmask.data(), coolunit, gasgr.data(),
        myisrf.data(), kappa_tot.data(), my_chemistry, my_rates, my_fields,
        idx_range, grain_temperatures, gas_grainsp_heatrate, grain_kappa,
        logTlininterp_buf, internal_dust_prop_buf
      );

      // Copy slice values back to grid

      for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
        if (itmask[ip1-1] != MASK_FALSE)  {
          if (my_chemistry->use_multiple_dust_temperatures == 0)  {
            dust_temp(ip1-1,j,k) = tdust[ip1-1];
          } else {
            if (my_chemistry->dust_species > 0)  {
              MgSiO3_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][ip1-1];
              AC_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::AC_dust][ip1-1];
            }
            if (my_chemistry->dust_species > 1)  {
              SiM_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][ip1-1];
              FeM_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][ip1-1];
              Mg2SiO4_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][ip1-1];
              Fe3O4_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][ip1-1];
              SiO2D_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][ip1-1];
              MgO_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][ip1-1];
              FeS_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][ip1-1];
              Al2O3_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][ip1-1];
            }
            if (my_chemistry->dust_species > 2)  {
              reforg_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][ip1-1];
              volorg_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][ip1-1];
              H2Oice_temp(ip1-1,j,k) = grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][ip1-1];
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

    // parallelize the k and j loops with OpenMP
    // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: !$omp parallel do schedule(runtime)
  //_// PORT: #endif
    //_// TODO_USE: OMP_PRAGMA("omp parallel")
    {
      //_// TODO: move relevant variable declarations to here to replace OMP private
      //_// TODO_USE: OMP_PRAGMA("omp for")
      for (int t = 0; t < idx_helper.outer_ind_size; t++) {
        // construct an index-range corresponding to "i-slice"
        const IndexRange idx_range = make_idx_range_(t, &idx_helper);
        const int k = idx_range.k;
        const int j = idx_range.j;

        if (imetal == 1)  {
          for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
            metal(ip1-1,j,k) = metal(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
            if (my_chemistry->multi_metals > 0)  {
              metal_loc(ip1-1,j,k) = metal_loc(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C13(ip1-1,j,k) = metal_C13(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C20(ip1-1,j,k) = metal_C20(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C25(ip1-1,j,k) = metal_C25(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C30(ip1-1,j,k) = metal_C30(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F13(ip1-1,j,k) = metal_F13(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F15(ip1-1,j,k) = metal_F15(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F50(ip1-1,j,k) = metal_F50(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F80(ip1-1,j,k) = metal_F80(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_P170(ip1-1,j,k)=metal_P170(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_P200(ip1-1,j,k)=metal_P200(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_Y19(ip1-1,j,k) = metal_Y19(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
            }
          }
        }
        if (my_chemistry->use_dust_density_field == 1)  {
          for (int ip1 = idx_range.i_start + 1; ip1<=(idx_range.i_end + 1); ip1++) {
            dust(ip1-1,j,k) = dust(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              // !            if (metal(ip1,j,k) .gt. 1.d-9 * d(ip1,j,k)) then
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(ip1-1,j,k)  = MgSiO3(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                AC(ip1-1,j,k)      = AC(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(ip1-1,j,k)     = SiM(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                FeM(ip1-1,j,k)     = FeM(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                Mg2SiO4(ip1-1,j,k) = Mg2SiO4(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                Fe3O4(ip1-1,j,k)   = Fe3O4(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                SiO2D(ip1-1,j,k)   = SiO2D(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                MgO(ip1-1,j,k)     = MgO(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                FeS(ip1-1,j,k)     = FeS(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                Al2O3(ip1-1,j,k)   = Al2O3(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(ip1-1,j,k)  = reforg(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                volorg(ip1-1,j,k)  = volorg(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
                H2Oice(ip1-1,j,k)  = H2Oice(ip1-1,j,k)*(gr_float)(std::pow(internalu.a_value,3) );
              }
              // !            endif
            }
          }
        }
      }
    }  // OMP_PRAGMA("omp parallel")

  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
