// See LICENSE file for license and copyright information

/// @file calc_tdust_3d.C
/// @brief Declares signature of calc_tdust_3d_g

// This file was initially generated automatically during conversion of the
// calc_tdust_3d_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "calc_tdust_3d.h"
#include "fortran_func_decls.h"
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

  // grain growth
  std::vector<double> sgSiM(my_fields->grid_dimension[0]);
  std::vector<double> sgFeM(my_fields->grid_dimension[0]);
  std::vector<double> sgMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> sgMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> sgFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> sgAC(my_fields->grid_dimension[0]);
  std::vector<double> sgSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> sgMgO(my_fields->grid_dimension[0]);
  std::vector<double> sgFeS(my_fields->grid_dimension[0]);
  std::vector<double> sgAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> sgreforg(my_fields->grid_dimension[0]);
  std::vector<double> sgvolorg(my_fields->grid_dimension[0]);
  std::vector<double> sgH2Oice(my_fields->grid_dimension[0]);
  std::vector<double> sgtot(my_fields->grid_dimension[0]);
  std::vector<double> alSiM(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFeM(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMg2SiO4(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMgSiO3(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFe3O4(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alAC(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alSiO2D(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alMgO(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alFeS(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alAl2O3(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alreforg(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alvolorg(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> alH2Oice(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> altot(my_rates->gr_N[2-1] * my_fields->grid_dimension[0]);
  std::vector<double> kpSiM(my_fields->grid_dimension[0]);
  std::vector<double> kpFeM(my_fields->grid_dimension[0]);
  std::vector<double> kpMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> kpMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> kpFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> kpAC(my_fields->grid_dimension[0]);
  std::vector<double> kpSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> kpMgO(my_fields->grid_dimension[0]);
  std::vector<double> kpFeS(my_fields->grid_dimension[0]);
  std::vector<double> kpAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> kpreforg(my_fields->grid_dimension[0]);
  std::vector<double> kpvolorg(my_fields->grid_dimension[0]);
  std::vector<double> kpH2Oice(my_fields->grid_dimension[0]);
  std::vector<double> kptot(my_fields->grid_dimension[0]);

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
  // -- removed line (previously just declared arg types) -- 
  std::vector<double> gasSiM(my_fields->grid_dimension[0]);
  std::vector<double> gasFeM(my_fields->grid_dimension[0]);
  std::vector<double> gasMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> gasMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> gasFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> gasAC(my_fields->grid_dimension[0]);
  std::vector<double> gasSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> gasMgO(my_fields->grid_dimension[0]);
  std::vector<double> gasFeS(my_fields->grid_dimension[0]);
  std::vector<double> gasAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> gasreforg(my_fields->grid_dimension[0]);
  std::vector<double> gasvolorg(my_fields->grid_dimension[0]);
  std::vector<double> gasH2Oice(my_fields->grid_dimension[0]);

  // Parameters

  const double mh_local_var = mh_grflt;

  // Locals

  double trad, zr, logtem0, logtem9, dlogtem, coolunit, dbase1, tbase1, xbase1;

  // Slice locals
 
  std::vector<long long> indixe(my_fields->grid_dimension[0]);
  std::vector<double> t1(my_fields->grid_dimension[0]);
  std::vector<double> t2(my_fields->grid_dimension[0]);
  std::vector<double> logtem(my_fields->grid_dimension[0]);
  std::vector<double> tdef(my_fields->grid_dimension[0]);
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
    //_// PORT: !$omp parallel do schedule(runtime) private(i, j, k)
    //_// PORT: #endif
    //_// TODO_USE: OMP_PRAGMA("omp parallel")
    {
      //_// TODO: move relevant variable declarations to here to replace OMP private

      //_// TODO_USE: OMP_PRAGMA("omp for")
      for (int t = 0; t < idx_helper.outer_ind_size; t++) {
        // construct an index-range corresponding to "i-slice"
        const IndexRange idx_range = make_idx_range_(t, &idx_helper);
        int k = idx_range.kp1; // use 1-based indexing (for now)
        int j = idx_range.jp1; // use 1-based indexing (for now)

        if (imetal == 1)  {
          for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
            metal(i-1,j-1,k-1) = metal(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
            if (my_chemistry->multi_metals > 0)  {
              metal_loc(i-1,j-1,k-1) = metal_loc(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C13(i-1,j-1,k-1) = metal_C13(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C20(i-1,j-1,k-1) = metal_C20(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C25(i-1,j-1,k-1) = metal_C25(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_C30(i-1,j-1,k-1) = metal_C30(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F13(i-1,j-1,k-1) = metal_F13(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F15(i-1,j-1,k-1) = metal_F15(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F50(i-1,j-1,k-1) = metal_F50(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_F80(i-1,j-1,k-1) = metal_F80(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_P170(i-1,j-1,k-1)=metal_P170(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_P200(i-1,j-1,k-1)=metal_P200(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              metal_Y19(i-1,j-1,k-1) = metal_Y19(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
            }
          }
        }
        if (my_chemistry->use_dust_density_field == 1)  {
          for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
            dust(i-1,j-1,k-1) = dust(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              // !            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1)  = MgSiO3(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                AC(i-1,j-1,k-1)      = AC(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(i-1,j-1,k-1)     = SiM(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                FeM(i-1,j-1,k-1)     = FeM(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                Fe3O4(i-1,j-1,k-1)   = Fe3O4(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                SiO2D(i-1,j-1,k-1)   = SiO2D(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                MgO(i-1,j-1,k-1)     = MgO(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                FeS(i-1,j-1,k-1)     = FeS(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                Al2O3(i-1,j-1,k-1)   = Al2O3(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1)  = reforg(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                volorg(i-1,j-1,k-1)  = volorg(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
                H2Oice(i-1,j-1,k-1)  = H2Oice(i-1,j-1,k-1)/(gr_float)(std::pow(internalu.a_value,3) );
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
  //_// PORT: !$omp&   i, j, k,
  //_// PORT: !$omp&   indixe,
  //_// PORT: !$omp&   t1, t2, logtem, tdef,
  //_// PORT: !$omp&   tgas, tdust, nh, gasgr, myisrf,
  //_// PORT: !$omp&   itmask )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    //_// TODO: move relevant variable declarations to here to replace OMP private
    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      IndexRange idx_range = make_idx_range_(t, &idx_helper);
      int k = idx_range.kp1; // use 1-based indexing (for now)
      int j = idx_range.jp1; // use 1-based indexing (for now)

      for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {

        // Set itmask to all true

        itmask[i-1] = MASK_TRUE;

      }

      // Iteration mask for metal-rich cells

      if (imetal == 1)  {
        for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
          if (metal(i-1,j-1,k-1) < 1.e-9 * d(i-1,j-1,k-1))  {
            itmask[i-1] = MASK_FALSE;
          }
        }
      }

      // Compute grain size increment

      if ( (my_chemistry->use_dust_density_field > 0)  &&  (my_chemistry->dust_species > 0) )  {

         FORTRAN_NAME(calc_grain_size_increment_1d)(
                  &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->grain_growth, itmask.data(),
                 &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &idx_range.jp1, &idx_range.kp1, &dom, d.data(),
                 SiM.data(), FeM.data(), Mg2SiO4.data(), MgSiO3.data(), Fe3O4.data(),
                 AC.data(), SiO2D.data(), MgO.data(), FeS.data(), Al2O3.data(),
                 reforg.data(), volorg.data(), H2Oice.data(),
                 metal.data(), metal_loc.data(),
                 metal_C13.data(), metal_C20.data(), metal_C25.data(), metal_C30.data(),
                 metal_F13.data(), metal_F15.data(), metal_F50.data(), metal_F80.data(),
                 metal_P170.data(), metal_P200.data(), metal_Y19.data(),
                 &my_rates->SN0_N,
                 my_rates->SN0_fSiM, my_rates->SN0_fFeM, my_rates->SN0_fMg2SiO4, my_rates->SN0_fMgSiO3,
                 my_rates->SN0_fFe3O4, my_rates->SN0_fAC, my_rates->SN0_fSiO2D, my_rates->SN0_fMgO,
                 my_rates->SN0_fFeS, my_rates->SN0_fAl2O3,
                 my_rates->SN0_freforg, my_rates->SN0_fvolorg, my_rates->SN0_fH2Oice,
                 my_rates->SN0_r0SiM, my_rates->SN0_r0FeM, my_rates->SN0_r0Mg2SiO4, my_rates->SN0_r0MgSiO3,
                 my_rates->SN0_r0Fe3O4, my_rates->SN0_r0AC, my_rates->SN0_r0SiO2D, my_rates->SN0_r0MgO,
                 my_rates->SN0_r0FeS, my_rates->SN0_r0Al2O3,
                 my_rates->SN0_r0reforg, my_rates->SN0_r0volorg, my_rates->SN0_r0H2Oice,
                 my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td,
                 my_rates->SN0_kpSiM, my_rates->SN0_kpFeM, my_rates->SN0_kpMg2SiO4, my_rates->SN0_kpMgSiO3,
                 my_rates->SN0_kpFe3O4, my_rates->SN0_kpAC, my_rates->SN0_kpSiO2D, my_rates->SN0_kpMgO,
                 my_rates->SN0_kpFeS, my_rates->SN0_kpAl2O3,
                 my_rates->SN0_kpreforg, my_rates->SN0_kpvolorg, my_rates->SN0_kpH2Oice,
                 sgSiM.data(), sgFeM.data(), sgMg2SiO4.data(), sgMgSiO3.data(), sgFe3O4.data(), sgAC.data(),
                 sgSiO2D.data(), sgMgO.data(), sgFeS.data(), sgAl2O3.data(),
                 sgreforg.data(), sgvolorg.data(), sgH2Oice.data(), sgtot.data(),
                 alSiM.data(), alFeM.data(), alMg2SiO4.data(), alMgSiO3.data(), alFe3O4.data(), alAC.data(),
                 alSiO2D.data(), alMgO.data(), alFeS.data(), alAl2O3.data(),
                 alreforg.data(), alvolorg.data(), alH2Oice.data(), altot.data()
              );

      }

      for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if(itmask[i-1] != MASK_FALSE)  {
          // Calculate metallicity

          if (imetal == 1)  {
            metallicity[i-1] = metal(i-1,j-1,k-1) / d(i-1,j-1,k-1) / my_chemistry->SolarMetalFractionByMass;
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
            dust2gas[i-1] = dust(i-1,j-1,k-1) / d(i-1,j-1,k-1);
          } else {
            dust2gas[i-1] = my_chemistry->local_dust_to_gas_ratio * metallicity[i-1];
          }

          // Compute interstellar radiation field

          if (my_chemistry->use_isrf_field > 0)  {
            myisrf[i-1] = isrf_habing(i-1,j-1,k-1);
          } else {
            myisrf[i-1] = my_chemistry->interstellar_radiation_field;
          }


          // Compute hydrogen number density

          nh[i-1] = HI(i-1,j-1,k-1) + HII(i-1,j-1,k-1);
          if (my_chemistry->primordial_chemistry > 1)  {
            nh[i-1] = nh[i-1] + H2I(i-1,j-1,k-1) + H2II(i-1,j-1,k-1);
          }

          // We have not converted to proper, so use urho and not dom

          nh[i-1] = nh[i-1] * internalu.urho / mh_local_var;

          // Compute log temperature and truncate if above/below table max/min

          tgas[i-1]   = gas_temp(i-1,j-1,k-1);
          logtem[i-1] = std::log(tgas[i-1]);
          logtem[i-1] = std::fmax(logtem[i-1], logtem0);
          logtem[i-1] = std::fmin(logtem[i-1], logtem9);

          // Compute index into the table and precompute parts of linear interp

          indixe[i-1] = std::fmin(my_chemistry->NumberOfTemperatureBins-1,
               std::fmax(1,(long long)((logtem[i-1]-logtem0)/dlogtem)+1));
          t1[i-1] = (logtem0 + (indixe[i-1] - 1)*dlogtem);
          t2[i-1] = (logtem0 + (indixe[i-1]    )*dlogtem);
          tdef[i-1] = (logtem[i-1] - t1[i-1]) / (t2[i-1] - t1[i-1]);

        }
      }
      // --- Compute dust temperature in a slice ---

       FORTRAN_NAME(calc_all_tdust_gasgr_1d_g)(&my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                    &my_chemistry->use_dust_density_field, &my_fields->grid_start[0], &my_fields->grid_end[0], &idx_range.jp1, &idx_range.kp1, &my_chemistry->local_dust_to_gas_ratio, &my_rates->gamma_isrf,
                    &trad, my_rates->gas_grain, indixe.data(), tdef.data(), tgas.data(), tdust.data(),
                    metallicity.data(), dust2gas.data(), nh.data(), gasgr_tdust.data(),
                    itmask.data(),
                    &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT,
                    my_rates->gr_Td, grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
                    grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust],
                    grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], my_rates->gas_grain2, &my_rates->gamma_isrf2,
                    &coolunit, gasgr.data(), myisrf.data(), sgSiM.data(), sgFeM.data(), sgMg2SiO4.data(),
                    sgMgSiO3.data(), sgFe3O4.data(), sgAC.data(), sgSiO2D.data(), sgMgO.data(), sgFeS.data(),
                    sgAl2O3.data(), sgreforg.data(), sgvolorg.data(), sgH2Oice.data(), sgtot.data(),
                    alSiM.data(), alFeM.data(), alMg2SiO4.data(), alMgSiO3.data(), alFe3O4.data(), alAC.data(),
                    alSiO2D.data(), alMgO.data(), alFeS.data(), alAl2O3.data(), alreforg.data(),
                    alvolorg.data(), alH2Oice.data(), altot.data(), kpSiM.data(), kpFeM.data(),
                    kpMg2SiO4.data(), kpMgSiO3.data(), kpFe3O4.data(), kpAC.data(), kpSiO2D.data(),
                    kpMgO.data(), kpFeS.data(), kpAl2O3.data(), kpreforg.data(), kpvolorg.data(),
                    kpH2Oice.data(), kptot.data(), gasSiM.data(), gasFeM.data(), gasMg2SiO4.data(),
                    gasMgSiO3.data(), gasFe3O4.data(), gasAC.data(), gasSiO2D.data(), gasMgO.data(),
                    gasFeS.data(), gasAl2O3.data(), gasreforg.data(), gasvolorg.data(),
                    gasH2Oice.data()
                  );


      // Copy slice values back to grid

      for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
        if (itmask[i-1] != MASK_FALSE)  {
          if (my_chemistry->use_multiple_dust_temperatures == 0)  {
            dust_temp(i-1,j-1,k-1) = tdust[i-1];
          } else {
            if (my_chemistry->dust_species > 0)  {
              MgSiO3_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust][i-1];
              AC_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::AC_dust][i-1];
            }
            if (my_chemistry->dust_species > 1)  {
              SiM_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::SiM_dust][i-1];
              FeM_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::FeM_dust][i-1];
              Mg2SiO4_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust][i-1];
              Fe3O4_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust][i-1];
              SiO2D_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust][i-1];
              MgO_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::MgO_dust][i-1];
              FeS_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::FeS_dust][i-1];
              Al2O3_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust][i-1];
            }
            if (my_chemistry->dust_species > 2)  {
              reforg_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust][i-1];
              volorg_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust][i-1];
              H2Oice_temp(i-1,j-1,k-1) = grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust][i-1];
            }
          }
        }
      }

    }

    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from 'proper'

  if (internalu.extfields_in_comoving == 1)  {

    // parallelize the k and j loops with OpenMP
    // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: !$omp parallel do schedule(runtime) private(i, j, k)
  //_// PORT: #endif
    //_// TODO_USE: OMP_PRAGMA("omp parallel")
    {
      //_// TODO: move relevant variable declarations to here to replace OMP private
      //_// TODO_USE: OMP_PRAGMA("omp for")
      for (int t = 0; t < idx_helper.outer_ind_size; t++) {
        // construct an index-range corresponding to "i-slice"
        const IndexRange idx_range = make_idx_range_(t, &idx_helper);
        int k = idx_range.kp1; // use 1-based indexing (for now)
        int j = idx_range.jp1; // use 1-based indexing (for now)

        if (imetal == 1)  {
          for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
            metal(i-1,j-1,k-1) = metal(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
            if (my_chemistry->multi_metals > 0)  {
              metal_loc(i-1,j-1,k-1) = metal_loc(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C13(i-1,j-1,k-1) = metal_C13(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C20(i-1,j-1,k-1) = metal_C20(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C25(i-1,j-1,k-1) = metal_C25(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_C30(i-1,j-1,k-1) = metal_C30(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F13(i-1,j-1,k-1) = metal_F13(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F15(i-1,j-1,k-1) = metal_F15(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F50(i-1,j-1,k-1) = metal_F50(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_F80(i-1,j-1,k-1) = metal_F80(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_P170(i-1,j-1,k-1)=metal_P170(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_P200(i-1,j-1,k-1)=metal_P200(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              metal_Y19(i-1,j-1,k-1) = metal_Y19(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
            }
          }
        }
        if (my_chemistry->use_dust_density_field == 1)  {
          for (int i = idx_range.i_start + 1; i<=(idx_range.i_end + 1); i++) {
            dust(i-1,j-1,k-1) = dust(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
            if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
              // !            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (my_chemistry->dust_species > 0)  {
                MgSiO3(i-1,j-1,k-1)  = MgSiO3(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                AC(i-1,j-1,k-1)      = AC(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 1)  {
                SiM(i-1,j-1,k-1)     = SiM(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                FeM(i-1,j-1,k-1)     = FeM(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                Mg2SiO4(i-1,j-1,k-1) = Mg2SiO4(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                Fe3O4(i-1,j-1,k-1)   = Fe3O4(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                SiO2D(i-1,j-1,k-1)   = SiO2D(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                MgO(i-1,j-1,k-1)     = MgO(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                FeS(i-1,j-1,k-1)     = FeS(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                Al2O3(i-1,j-1,k-1)   = Al2O3(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
              }
              if (my_chemistry->dust_species > 2)  {
                reforg(i-1,j-1,k-1)  = reforg(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                volorg(i-1,j-1,k-1)  = volorg(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
                H2Oice(i-1,j-1,k-1)  = H2Oice(i-1,j-1,k-1)*(gr_float)(std::pow(internalu.a_value,3) );
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
