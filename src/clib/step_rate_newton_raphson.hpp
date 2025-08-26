// See LICENSE file for license and copyright information

/// @file step_rate_newton_raphson.hpp
/// @brief Defines the step_rate_newton_raphson function

// This file was initially generated automatically during conversion of the
// step_rate_newton_raphson function from FORTRAN to C++

#ifndef STEP_RATE_NEWTON_RAPHSON_HPP
#define STEP_RATE_NEWTON_RAPHSON_HPP

#include <vector>

#include "grackle.h"             // gr_float
#include "fortran_func_decls.h"  // gr_mask_int
//#include "fortran_func_wrappers.hpp" // grackle::impl::fortran_wrapper::gaussj_g
#include "gaussj_g.hpp"
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "utils-cpp.hpp"

#include "utils-field.hpp"
#include "time_deriv_0d.hpp"

namespace grackle::impl {

int sanity_check_() {
  GRIMPL_REQUIRE(((1+SpLUT::e)== 1), "index of e is %d", SpLUT::e);
  GRIMPL_REQUIRE(((1+SpLUT::HI)== 2), "index of HI is %d", SpLUT::HI);
  GRIMPL_REQUIRE(((1+SpLUT::HII)== 3), "index of HII is %d", SpLUT::HII);
  GRIMPL_REQUIRE(((1+SpLUT::HeI)== 4), "index of HeI is %d", SpLUT::HeI);
  GRIMPL_REQUIRE(((1+SpLUT::HeII)== 5), "index of HeII is %d", SpLUT::HeII);
  GRIMPL_REQUIRE(((1+SpLUT::HeIII)== 6), "index of HeIII is %d", SpLUT::HeIII);
  GRIMPL_REQUIRE(((1+SpLUT::HM)== 7), "index of HM is %d", SpLUT::HM);
  GRIMPL_REQUIRE(((1+SpLUT::H2I)== 8), "index of H2I is %d", SpLUT::H2I);
  GRIMPL_REQUIRE(((1+SpLUT::H2II)== 9), "index of H2II is %d", SpLUT::H2II);
  GRIMPL_REQUIRE(((1+SpLUT::DI)==10), "index of DI is %d", SpLUT::DI);
  GRIMPL_REQUIRE(((1+SpLUT::DII)==11), "index of DII is %d", SpLUT::DII);
  GRIMPL_REQUIRE(((1+SpLUT::HDI)==12), "index of HDI is %d", SpLUT::HDI);
  GRIMPL_REQUIRE(((1+SpLUT::DM)==13), "index of DM is %d", SpLUT::DM);
  GRIMPL_REQUIRE(((1+SpLUT::HDII)==14), "index of HDII is %d", SpLUT::HDII);
  GRIMPL_REQUIRE(((1+SpLUT::HeHII)==15), "index of HeHII is %d", SpLUT::HeHII);
  GRIMPL_REQUIRE(((1+SpLUT::CI)==16), "index of CI is %d", SpLUT::CI);
  GRIMPL_REQUIRE(((1+SpLUT::CII)==17), "index of CII is %d", SpLUT::CII);
  GRIMPL_REQUIRE(((1+SpLUT::CO)==18), "index of CO is %d", SpLUT::CO);
  GRIMPL_REQUIRE(((1+SpLUT::CO2)==19), "index of CO2 is %d", SpLUT::CO2);
  GRIMPL_REQUIRE(((1+SpLUT::OI)==20), "index of OI is %d", SpLUT::OI);
  GRIMPL_REQUIRE(((1+SpLUT::OH)==21), "index of OH is %d", SpLUT::OH);
  GRIMPL_REQUIRE(((1+SpLUT::H2O)==22), "index of H2O is %d", SpLUT::H2O);
  GRIMPL_REQUIRE(((1+SpLUT::O2)==23), "index of O2 is %d", SpLUT::O2);
  GRIMPL_REQUIRE(((1+SpLUT::SiI)==24), "index of SiI is %d", SpLUT::SiI);
  GRIMPL_REQUIRE(((1+SpLUT::SiOI)==25), "index of SiOI is %d", SpLUT::SiOI);
  GRIMPL_REQUIRE(((1+SpLUT::SiO2I)==26), "index of SiO2I is %d", SpLUT::SiO2I);
  GRIMPL_REQUIRE(((1+SpLUT::CH)==27), "index of CH is %d", SpLUT::CH);
  GRIMPL_REQUIRE(((1+SpLUT::CH2)==28), "index of CH2 is %d", SpLUT::CH2);
  GRIMPL_REQUIRE(((1+SpLUT::COII)==29), "index of COII is %d", SpLUT::COII);
  GRIMPL_REQUIRE(((1+SpLUT::OII)==30), "index of OII is %d", SpLUT::OII);
  GRIMPL_REQUIRE(((1+SpLUT::OHII)==31), "index of OHII is %d", SpLUT::OHII);
  GRIMPL_REQUIRE(((1+SpLUT::H2OII)==32), "index of H2OII is %d", SpLUT::H2OII);
  GRIMPL_REQUIRE(((1+SpLUT::H3OII)==33), "index of H3OII is %d", SpLUT::H3OII);
  GRIMPL_REQUIRE(((1+SpLUT::O2II)==34), "index of O2II is %d", SpLUT::O2II);
  GRIMPL_REQUIRE(((1+SpLUT::Mg)==35), "index of Mg is %d", SpLUT::Mg);
  GRIMPL_REQUIRE(((1+SpLUT::Al)==36), "index of Al is %d", SpLUT::Al);
  GRIMPL_REQUIRE(((1+SpLUT::S)==37), "index of S is %d", SpLUT::S);
  GRIMPL_REQUIRE(((1+SpLUT::Fe)==38), "index of Fe is %d", SpLUT::Fe);
  GRIMPL_REQUIRE(((1+SpLUT::MgSiO3_dust)==39), "index of MgSiO3 is %d", SpLUT::MgSiO3_dust);
  GRIMPL_REQUIRE(((1+SpLUT::AC_dust)==40), "index of AC is %d", SpLUT::AC_dust);
  GRIMPL_REQUIRE(((1+SpLUT::SiM_dust)==41), "index of SiM is %d", SpLUT::SiM_dust);
  GRIMPL_REQUIRE(((1+SpLUT::FeM_dust)==42), "index of FeM is %d", SpLUT::FeM_dust);
  GRIMPL_REQUIRE(((1+SpLUT::Mg2SiO4_dust)==43), "index of Mg2SiO4 is %d", SpLUT::Mg2SiO4_dust);
  GRIMPL_REQUIRE(((1+SpLUT::Fe3O4_dust)==44), "index of Fe3O4 is %d", SpLUT::Fe3O4_dust);
  GRIMPL_REQUIRE(((1+SpLUT::SiO2_dust)==45), "index of SiO2_dust is %d", SpLUT::SiO2_dust);
  GRIMPL_REQUIRE(((1+SpLUT::MgO_dust)==46), "index of MgO is %d", SpLUT::MgO_dust);
  GRIMPL_REQUIRE(((1+SpLUT::FeS_dust)==47), "index of FeS is %d", SpLUT::FeS_dust);
  GRIMPL_REQUIRE(((1+SpLUT::Al2O3_dust)==48), "index of Al2O3 is %d", SpLUT::Al2O3_dust);
  GRIMPL_REQUIRE(((1+SpLUT::ref_org_dust)==49), "index of reforg is %d", SpLUT::ref_org_dust);
  GRIMPL_REQUIRE(((1+SpLUT::vol_org_dust)==50), "index of volorg is %d", SpLUT::vol_org_dust);
  GRIMPL_REQUIRE(((1+SpLUT::H2O_ice_dust)==51), "index of H2Oice is %d", SpLUT::H2O_ice_dust);

  //GRIMPL_REQUIRE(false, "early abort");

  return 1;
}

/// acts as a temporary adaptor layer for time_deriv_0d::derivatives:
///
/// @note
/// Some additional context
/// - during transcription, time_deriv_0d::derivatives (formerly called
///   lookup_cool_rates0d) needed to be adapted to make use of data structures
///   already adopted in other parts of the transcribed code
///   - (recall that these data structures were introduced to remove HUNDREDS
///     of function arguments and HUNDREDS of local variables)
/// - since the newton-raphson solver already does a bunch of "translation"
///   between data representations, it should directly map into these preferred
///   data structures
/// - as a short-term (and extremely temporary) solution, we this function
///   to copy between formats & wrap the time_deriv_0d::derivatives function.
///   This allows to finish polishing time_deriv_0d::derivatives without
///   refactoring the newton-raphson scheme)
inline void wrapped_calc_derivatives(
  double dt,
  const double* dsp,
  double* dspdot,
  grackle::impl::time_deriv_0d::ContextPack& pack,
  std::vector<gr_float>& rhosp_grflt,         // <- preallocated temp buffer
  grackle::impl::SpeciesCollection& rhosp_dot // <- preallocated temp buffer
) {

  gr_float local_eint[1];
  double eint_dot_specific[1];

  // 1. copy values out of dsp into rhosp_grflt and local_eint
  //    - at some point, we need to add some logic to deal with edge cases
  //      that can arise when `sizeof(gr_float) < sizeof(double)
  for (int sp_idx = 0; sp_idx < SpLUT::NUM_ENTRIES; sp_idx++) {
    rhosp_grflt[sp_idx] = static_cast<gr_float>(dsp[sp_idx]);
  }
  local_eint[0] = static_cast<gr_float>(dsp[SpLUT::NUM_ENTRIES]);

  // 2. call the wrapped function
  grackle::impl::time_deriv_0d::derivatives(
    dt, rhosp_grflt.data(), rhosp_dot, local_eint, eint_dot_specific, pack
  );

  // 3. copy the computed derivatives into dspdot
  for (int sp_idx = 0; sp_idx < SpLUT::NUM_ENTRIES; sp_idx++) {
    dspdot[sp_idx] = rhosp_dot.data[sp_idx][0];
  }
  dspdot[SpLUT::NUM_ENTRIES] = eint_dot_specific[0];
}

/// An alternative to step_rate_g for evolving the species rate equations that
/// employs the Newton-Raphson scheme rather than Gauss-Seidel
///
/// Whereas the Gauss-Seidel scheme always operator-splits evolution of the
/// energy equation, this function internally evolves the energy equation
///
/// @note
/// The values in the various buffers (p2d, tgas, tdust, metallicity, dust2gas,
/// rhoH, mmw, h2dust, edot, grain_temperatures, logTlininterp_buf,
/// cool1dmulti_buf, coolingheating_buf, chemheatrates_buf), hold undefined
/// values after this function call wherever `itmask_nr` indicates that this
/// function runs.
/// - this has **ALWAYS** been the case. Historically these buffers have been
///   reused when computing finite differences
inline void step_rate_newton_raphson(
  int imetal, IndexRange idx_range, int iter, double dom, double chunit,
  double dx_cgs, double c_ljeans, double* dtit, double* p2d, double* tgas,
  double* tdust, double* metallicity, double* dust2gas, double* rhoH,
  double* mmw, double* h2dust, double* edot, gr_mask_type anydust,
  gr_mask_type* itmask_nr, gr_mask_type* itmask_metal, int* imp_eng,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates,
  InternalGrUnits internalu,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf,
  grackle::impl::CoolHeatScratchBuf coolingheating_buf,
  grackle::impl::ChemHeatingRates chemheatrates_buf
)
{
  // shorten `grackle::impl::time_deriv_0d` to `t_deriv` within this function
  namespace t_deriv = ::grackle::impl::time_deriv_0d;
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  // the following is a quick sanity check
  // -> we can remove this once we finish transcribing lookup_cool_rate0d
  sanity_check_();

  // Density, energy and velocity fields fields

  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // ierror local variable
  //   - this variable is only used internally by this subroutine for
  //     determining control-flow
  //   - it does NOT report whether or not this function succeeded (maybe
  //     we should change the name?)
  int ierror;
  // Local variable
  int itr, itr_time;
  int nsp, isp, jsp, id;
  double dspj, err, err_max;
  // the following specifies the historical 1-based index that we would use to
  // hold energy
  const int i_eng = SpLUT::NUM_ENTRIES + 1;
  // There may be an argument for allocating the following at a higher
  // level function, but we will leave that for after transcription
  std::vector<double> dsp(i_eng);
  std::vector<double> dsp0(i_eng);
  std::vector<double> dsp1(i_eng);
  std::vector<double> dspdot(i_eng);
  std::vector<double> dspdot1(i_eng);
  std::vector<double> ddsp(i_eng);
  std::vector<double> der_data_(i_eng * i_eng);
  grackle::impl::View<double**> der(der_data_.data(), i_eng, i_eng);

  // (In the future, we may want to reconsider when/how we allocate
  // the following 3 variables)
  std::vector<int> idsp;
  std::vector<double> mtrx_data_;
  grackle::impl::View<double**> mtrx;
  std::vector<double> vec;
  // Another parameter
  const double eps = 1.e-4;

  // allocate some scratch buffers. In the future, when we refactor to support
  // calculations of the time derivatives of multiple values at once, we should
  // forward the buffers passed into this routine as arguments rather than
  // allocating separate buffers
  t_deriv::MainScratchBuf main_scratch_buf = t_deriv::new_MainScratchBuf();

  // collect args that are forwarded to the time-derivative calculation and are
  // effectively frozen between various calls
  t_deriv::FrozenSimpleArgs frozen_tderiv_args = {
    imetal, iter, dom, chunit, dx_cgs, c_ljeans, anydust, my_chemistry,
    my_rates, my_uvb_rates, internalu
  };

  // partially initialize the struct we will use for the time derivative calc
  t_deriv::ContextPack pack = t_deriv::new_ContextPack(
    frozen_tderiv_args, main_scratch_buf
  );

  std::vector<gr_float> rhosp_grflt(SpLUT::NUM_ENTRIES);
  grackle::impl::SpeciesCollection rhosp_dot =
    grackle::impl::new_SpeciesCollection(1);

  // the following check was inspired by a compiler warning indicating that
  // nsp won't be initialized if this condition isn't met
  GRIMPL_REQUIRE((my_chemistry->primordial_chemistry > 0),
    "this function can't support primordial_chemistry == %d",
    my_chemistry->primordial_chemistry
  );

  // The following was extracted from another subroutine
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    const int j = idx_range.j;
    const int k = idx_range.k;

    // remap the 3D index to a 1D index
    const int field_idx1d = layoutleft_3D_index_to_1D_(
        my_fields->grid_dimension, i, j, k
    );

    if (itmask_nr[i] != MASK_FALSE)  {

      // If density and temperature are low, update gas energy explicitly

      if (my_chemistry->with_radiative_cooling == 1)  {
        if (imp_eng[i] == 0)  {
          e(i,j,k)  = e(i,j,k) +
               (gr_float)(edot[i]/d(i,j,k)*dtit[i] );
        }
      }

      // initialize arrays
      nsp = 6; // (my_chemistry->primordial_chemistry >= 1)
      if (my_chemistry->primordial_chemistry > 1) { nsp = nsp + 3; }
      if (my_chemistry->primordial_chemistry > 2) { nsp = nsp + 3; }
      if (my_chemistry->primordial_chemistry > 3) { nsp = nsp + 3; }
      if (itmask_metal[i] != MASK_FALSE)  {
        if (my_chemistry->metal_chemistry == 1)  {
          nsp = nsp + 19;
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0) { nsp = nsp + 1; }
            if (my_chemistry->dust_species > 1) { nsp = nsp + 3; }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0) { nsp = nsp + 2; }
          if (my_chemistry->dust_species > 1) { nsp = nsp + 8; }
          if (my_chemistry->dust_species > 2) { nsp = nsp + 3; }
        }
      }
      nsp = nsp + imp_eng[i];
      idsp.reserve(nsp);
      mtrx_data_.reserve(nsp * nsp);
      mtrx = grackle::impl::View<double**>(mtrx_data_.data(), nsp, nsp);
      vec.reserve(nsp);

      // copy values into dsp from my_fields
      // -> in the future, we will be able write the next ~80 lines as
      //    a concise 4-line for-loop when we start to use the
      //    SpeciesLUTFieldAdaptor type

      if ( my_chemistry->primordial_chemistry > 0 )  {
        dsp[SpLUT::e] = my_fields->e_density[field_idx1d];
        dsp[SpLUT::HI] = my_fields->HI_density[field_idx1d];
        dsp[SpLUT::HII] = my_fields->HII_density[field_idx1d];
        dsp[SpLUT::HeI] = my_fields->HeI_density[field_idx1d];
        dsp[SpLUT::HeII] = my_fields->HeII_density[field_idx1d];
        dsp[SpLUT::HeIII] = my_fields->HeIII_density[field_idx1d];
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        dsp[SpLUT::HM] = my_fields->HM_density[field_idx1d];
        dsp[SpLUT::H2I] = my_fields->H2I_density[field_idx1d];
        dsp[SpLUT::H2II] = my_fields->H2II_density[field_idx1d];
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        dsp[SpLUT::DI] = my_fields->DI_density[field_idx1d];
        dsp[SpLUT::DII] = my_fields->DII_density[field_idx1d];
        dsp[SpLUT::HDI] = my_fields->HDI_density[field_idx1d];
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        dsp[SpLUT::DM] = my_fields->DM_density[field_idx1d];
        dsp[SpLUT::HDII] = my_fields->HDII_density[field_idx1d];
        dsp[SpLUT::HeHII] = my_fields->HeHII_density[field_idx1d];
      }
      if ( itmask_metal[i] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          dsp[SpLUT::CI] = my_fields->CI_density[field_idx1d];
          dsp[SpLUT::CII] = my_fields->CII_density[field_idx1d];
          dsp[SpLUT::CO] = my_fields->CO_density[field_idx1d];
          dsp[SpLUT::CO2] = my_fields->CO2_density[field_idx1d];
          dsp[SpLUT::OI] = my_fields->OI_density[field_idx1d];
          dsp[SpLUT::OH] = my_fields->OH_density[field_idx1d];
          dsp[SpLUT::H2O] = my_fields->H2O_density[field_idx1d];
          dsp[SpLUT::O2] = my_fields->O2_density[field_idx1d];
          dsp[SpLUT::SiI] = my_fields->SiI_density[field_idx1d];
          dsp[SpLUT::SiOI] = my_fields->SiOI_density[field_idx1d];
          dsp[SpLUT::SiO2I] = my_fields->SiO2I_density[field_idx1d];
          dsp[SpLUT::CH] = my_fields->CH_density[field_idx1d];
          dsp[SpLUT::CH2] = my_fields->CH2_density[field_idx1d];
          dsp[SpLUT::COII] = my_fields->COII_density[field_idx1d];
          dsp[SpLUT::OII] = my_fields->OII_density[field_idx1d];
          dsp[SpLUT::OHII] = my_fields->OHII_density[field_idx1d];
          dsp[SpLUT::H2OII] = my_fields->H2OII_density[field_idx1d];
          dsp[SpLUT::H3OII] = my_fields->H3OII_density[field_idx1d];
          dsp[SpLUT::O2II] = my_fields->O2II_density[field_idx1d];
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) ) {
            if (my_chemistry->dust_species > 0)  {
              dsp[SpLUT::Mg] = my_fields->Mg_density[field_idx1d];
            }
            if (my_chemistry->dust_species > 1)  {
              dsp[SpLUT::Al] = my_fields->Al_density[field_idx1d];
              dsp[SpLUT::S] = my_fields->S_density[field_idx1d];
              dsp[SpLUT::Fe] = my_fields->Fe_density[field_idx1d];
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) ) {
          if (my_chemistry->dust_species > 0)  {
            dsp[SpLUT::MgSiO3_dust] = my_fields->MgSiO3_dust_density[field_idx1d];
            dsp[SpLUT::AC_dust] = my_fields->AC_dust_density[field_idx1d];
          }
          if (my_chemistry->dust_species > 1)  {
            dsp[SpLUT::SiM_dust] = my_fields->SiM_dust_density[field_idx1d];
            dsp[SpLUT::FeM_dust] = my_fields->FeM_dust_density[field_idx1d];
            dsp[SpLUT::Mg2SiO4_dust] = my_fields->Mg2SiO4_dust_density[field_idx1d];
            dsp[SpLUT::Fe3O4_dust] = my_fields->Fe3O4_dust_density[field_idx1d];
            dsp[SpLUT::SiO2_dust] = my_fields->SiO2_dust_density[field_idx1d];
            dsp[SpLUT::MgO_dust] = my_fields->MgO_dust_density[field_idx1d];
            dsp[SpLUT::FeS_dust] = my_fields->FeS_dust_density[field_idx1d];
            dsp[SpLUT::Al2O3_dust] = my_fields->Al2O3_dust_density[field_idx1d];
          }
          if (my_chemistry->dust_species > 2)  {
            dsp[SpLUT::ref_org_dust] = my_fields->ref_org_dust_density[field_idx1d];
            dsp[SpLUT::vol_org_dust] = my_fields->vol_org_dust_density[field_idx1d];
            dsp[SpLUT::H2O_ice_dust] = my_fields->H2O_ice_dust_density[field_idx1d];
          }
        }
      }
      dsp[i_eng-1] = e(i,j,k);

      // here, we fill in the idsp array
      // -> the idsp array is used for mapping between the compressed vector
      //    form (that only has nsp elements) and dsp (that always has
      //    `SpLUT::NUM_ENTRIES + 1` entries)
      // -> in the future, we should construct this array first (before doing
      //    anything else):
      //    -> it gives us the value of nsp for free!
      //    -> using this array with SpeciesLUTFieldAdaptor lets us cut out
      //       (duplicated) logic when we copy values from my_fields to dsp and
      //       the reverse (we can cut ~70 lines in each place)
      id = 0;
      if ( my_chemistry->primordial_chemistry > 0 )  {
        idsp[id++] = SpLUT::e;
        idsp[id++] = SpLUT::HI;
        idsp[id++] = SpLUT::HII;
        idsp[id++] = SpLUT::HeI;
        idsp[id++] = SpLUT::HeII;
        idsp[id++] = SpLUT::HeIII;
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        idsp[id++] = SpLUT::HM;
        idsp[id++] = SpLUT::H2I;
        idsp[id++] = SpLUT::H2II;
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        idsp[id++] = SpLUT::DI;
        idsp[id++] = SpLUT::DII;
        idsp[id++] = SpLUT::HDI;
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        idsp[id++] = SpLUT::DM;
        idsp[id++] = SpLUT::HDII;
        idsp[id++] = SpLUT::HeHII;
      }
      if ( itmask_metal[i] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          idsp[id++] = SpLUT::CI;
          idsp[id++] = SpLUT::CII;
          idsp[id++] = SpLUT::CO;
          idsp[id++] = SpLUT::CO2;
          idsp[id++] = SpLUT::OI;
          idsp[id++] = SpLUT::OH;
          idsp[id++] = SpLUT::H2O;
          idsp[id++] = SpLUT::O2;
          idsp[id++] = SpLUT::SiI;
          idsp[id++] = SpLUT::SiOI;
          idsp[id++] = SpLUT::SiO2I;
          idsp[id++] = SpLUT::CH;
          idsp[id++] = SpLUT::CH2;
          idsp[id++] = SpLUT::COII;
          idsp[id++] = SpLUT::OII;
          idsp[id++] = SpLUT::OHII;
          idsp[id++] = SpLUT::H2OII;
          idsp[id++] = SpLUT::H3OII;
          idsp[id++] = SpLUT::O2II;
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
            if (my_chemistry->dust_species > 0)  {
              idsp[id++] = SpLUT::Mg;
            }
            if (my_chemistry->dust_species > 1)  {
              idsp[id++] = SpLUT::Al;
              idsp[id++] = SpLUT::S;
              idsp[id++] = SpLUT::Fe;
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
          if (my_chemistry->dust_species > 0)  {
            idsp[id++] = SpLUT::MgSiO3_dust;
            idsp[id++] = SpLUT::AC_dust;
          }
          if (my_chemistry->dust_species > 1)  {
            idsp[id++] = SpLUT::SiM_dust;
            idsp[id++] = SpLUT::FeM_dust;
            idsp[id++] = SpLUT::Mg2SiO4_dust;
            idsp[id++] = SpLUT::Fe3O4_dust;
            idsp[id++] = SpLUT::SiO2_dust;
            idsp[id++] = SpLUT::MgO_dust;
            idsp[id++] = SpLUT::FeS_dust;
            idsp[id++] = SpLUT::Al2O3_dust;
          }
          if (my_chemistry->dust_species > 2)  {
            idsp[id++] = SpLUT::ref_org_dust;
            idsp[id++] = SpLUT::vol_org_dust;
            idsp[id++] = SpLUT::H2O_ice_dust;
          }
        }
      }

      if ( imp_eng[i] ==1 )  {
        idsp[id++] = i_eng-1;
      }

      GRIMPL_REQUIRE((id == nsp), "SANITY CHECK FAILED");

      // setup the dt_FIXME variable (this obviously needs some attention)
      // -> as in solve_rate_cool_g, dt represents the total timestep that we
      //    are evolving the system of equations over
      // -> here, we are saying defining dt as a reference to the current
      //    subcycle timestep at the current zone (we define it as a
      //    reference because it needs to match the value of dtit[i] if we
      //    choose to modify it)
      // -> This sure seems inconsistent!
      // -> OK, in more detail:
      //    -> the only purpose of this variable is to be passed into the
      //       time-derivative calculation.
      //    -> within the time derivative calculation, the variable is ONLY
      //       used as an argument forwarded to the `dt` argument of
      //       lookup_cool_rates1d_g
      //    -> at present, the other callsite of `lookup_cool_rates1d_g`
      //       passes in the total timestep to dt
      //    -> for now, we're just going to maintain compatability but this
      //       NEEDS to be addressed (see the C++ docstring of
      //       `lookup_cool_rates1d_g` for more discussion). We actually have
      //       been doing things correctly here and incorrectly elsewhere
      double& dt_FIXME = dtit[i];

      // configure pack for use during the current index
      t_deriv::configure_ContextPack(
        &pack, my_fields, field_idx1d, itmask_metal[i], imp_eng[i]
      );

      // copy values at the current index from the external scratch buffers
      // passed into the current function as args into the corresponding
      // buffers managed by pack
      // -> see the called function's docstring for an explanation as to why
      //    this is not currently combined with the prior function call
      // -> in a lot of cases, the need to copy values is dictated by the value
      //    of imp_eng
      // -> in the future, we should really refactor so that this isn't
      //    necessary (arguably, the only thing we should really consider
      //    copying is the value of cool1dmulti_buf.tgasold -- and that doesn't
      //    currently get used)
      t_deriv::scratchbufs_copy_into_pack(
        i, &pack, p2d, tgas, tdust, metallicity, dust2gas, rhoH, mmw, h2dust,
        edot, grain_temperatures, logTlininterp_buf, cool1dmulti_buf,
        coolingheating_buf, chemheatrates_buf
      );

      // Save arrays at ttot(ip1)

      std::memcpy(dsp0.data(), dsp.data(), sizeof(double)*i_eng);
      std::memset(ddsp.data(), 0, sizeof(double)*i_eng);

      // Search for the timestep for which chemistry converges

      ierror=1;
      itr_time=0;
      while (ierror==1) {

        // If not converge, restore arrays at ttot(ip1)

        std::memcpy(dsp.data(), dsp0.data(), sizeof(double)*i_eng);
        std::memset(ddsp.data(), 0, sizeof(double)*i_eng);

        // Iteration to solve ODEs

        err_max=1.e2;
        itr=0;
        while (err_max>1.e-8) {
          if(itr>=20)  {
            ierror = 1;
            goto label_9996;
          }

          // calc the time derivatives
          wrapped_calc_derivatives(
            dt_FIXME, dsp.data(), dspdot.data(), pack, rhosp_grflt, rhosp_dot
          );

          for (jsp = 1; jsp<=(nsp); jsp++) {
            dspj = eps * dsp[idsp[jsp-1]];
            for (isp = 1; isp<=(nsp); isp++) {
              if(isp == jsp)  {
                dsp1[idsp[isp-1]] = dsp[idsp[isp-1]] + dspj;
              } else {
                dsp1[idsp[isp-1]] = dsp[idsp[isp-1]];
              }
            }

            wrapped_calc_derivatives(
              dt_FIXME, dsp1.data(), dspdot1.data(), pack, rhosp_grflt,
              rhosp_dot
            );

            for (isp = 1; isp<=(nsp); isp++) {
              if ( (dsp[idsp[isp-1]]==0.e0)
               &&  (dspdot1[idsp[isp-1]]
               ==  dspdot[idsp[isp-1]]) )  {
                der(idsp[isp-1],idsp[jsp-1]) = 0.e0;
              } else {
                der(idsp[isp-1],idsp[jsp-1]) =
                   (dspdot1[idsp[isp-1]]
                   - dspdot[idsp[isp-1]]) / dspj;
              }
            }

          }

          for (isp = 1; isp<=(nsp); isp++) {
            for (jsp = 1; jsp<=(nsp); jsp++) {
              if(isp == jsp)  {
                mtrx(isp-1,jsp-1) = 1.e0 - dtit[i]
                   * der(idsp[isp-1],idsp[jsp-1]);
              } else {
                mtrx(isp-1,jsp-1) =      - dtit[i]
                   * der(idsp[isp-1],idsp[jsp-1]);
              }
            }
          }

          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = dspdot[idsp[isp-1]] * dtit[i]
                   - ddsp[idsp[isp-1]];
          }

          // to get more accuracy
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]/d(i,j,k);
          }

          ierror = grackle::impl::gaussj_g(nsp, mtrx.data(), vec.data());
          if(ierror == 1)  {
            goto label_9998;
          }

          // multiply with density again
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]*d(i,j,k);
          }

          for (isp = 1; isp<=(nsp); isp++) {
            ddsp[idsp[isp-1]] = ddsp[idsp[isp-1]] + vec[isp-1];
            dsp[idsp[isp-1]]  = dsp[idsp[isp-1]]  + vec[isp-1];
          }

          if (imp_eng[i] == 1)  {
            if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
              for (isp = 1; isp<=(nsp); isp++) {
                if ( (dsp[idsp[isp-1]] != dsp[idsp[isp-1]])
                 ||  (dsp[idsp[isp-1]] <= 0.) )  {
                  ierror = 1;
                  goto label_9997;
                }
              }
            }
          }

          err_max = 0.e0;
          for (isp = 1; isp<=(nsp); isp++) {
            if(dsp[idsp[isp-1]] > tiny8)  {
              err = grackle::impl::dabs(vec[isp-1] / dsp[idsp[isp-1]]);
            } else {
              err = 0.e0;
            }
            if(err > err_max)  {
              err_max = err;
            }
          }

          itr=itr+1;
        }

label_9998:
label_9997:
label_9996:


        // Check if the fractions are valid after an iteration

        if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
          for (isp = 1; isp<=(nsp); isp++) {
            if ( (dsp[idsp[isp-1]] != dsp[idsp[isp-1]])
             ||  (dsp[idsp[isp-1]] <= 0.) )  {
              ierror = 1;
            }
          }
        }
        if(ierror == 1)  {
          dtit[i] = 0.5e0*dtit[i];
        }

        itr_time=itr_time+1;
      }


      // overwrite the scratch-buffer values
      // -> It makes no logical sense for us to do this, but we do it for the
      //    sake of backwards compatability
      // -> we should totally delete this function (we may already be able to
      //    do so)
      t_deriv::scratchbufs_copy_from_pack(
        i, &pack, p2d, tgas, tdust, metallicity, dust2gas, rhoH, mmw, h2dust,
        edot, grain_temperatures, logTlininterp_buf, cool1dmulti_buf,
        coolingheating_buf, chemheatrates_buf
      );

      // overwrite the fields tracked by my_fields with the evolved values
      // -> in the future, we will be able write the next ~80 lines as
      //    a concise 4-line for-loop when we start to use the
      //    SpeciesLUTFieldAdaptor type

      if ( my_chemistry->primordial_chemistry > 0 )  {
        my_fields->e_density[field_idx1d]       = dsp[SpLUT::e];
        my_fields->HI_density[field_idx1d]      = dsp[SpLUT::HI];
        my_fields->HII_density[field_idx1d]     = dsp[SpLUT::HII];
        my_fields->HeI_density[field_idx1d]     = dsp[SpLUT::HeI];
        my_fields->HeII_density[field_idx1d]    = dsp[SpLUT::HeII];
        my_fields->HeIII_density[field_idx1d]   = dsp[SpLUT::HeIII];
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        my_fields->HM_density[field_idx1d]      = dsp[SpLUT::HM];
        my_fields->H2I_density[field_idx1d]     = dsp[SpLUT::H2I];
        my_fields->H2II_density[field_idx1d]    = dsp[SpLUT::H2II];
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        my_fields->DI_density[field_idx1d]      = dsp[SpLUT::DI];
        my_fields->DII_density[field_idx1d]     = dsp[SpLUT::DII];
        my_fields->HDI_density[field_idx1d]     = dsp[SpLUT::HDI];
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        my_fields->DM_density[field_idx1d]      = dsp[SpLUT::DM];
        my_fields->HDII_density[field_idx1d]    = dsp[SpLUT::HDII];
        my_fields->HeHII_density[field_idx1d]   = dsp[SpLUT::HeHII];
      }
      if ( itmask_metal[i] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          my_fields->CI_density[field_idx1d]      = dsp[SpLUT::CI];
          my_fields->CII_density[field_idx1d]     = dsp[SpLUT::CII];
          my_fields->CO_density[field_idx1d]      = dsp[SpLUT::CO];
          my_fields->CO2_density[field_idx1d]     = dsp[SpLUT::CO2];
          my_fields->OI_density[field_idx1d]      = dsp[SpLUT::OI];
          my_fields->OH_density[field_idx1d]      = dsp[SpLUT::OH];
          my_fields->H2O_density[field_idx1d]     = dsp[SpLUT::H2O];
          my_fields->O2_density[field_idx1d]      = dsp[SpLUT::O2];
          my_fields->SiI_density[field_idx1d]     = dsp[SpLUT::SiI];
          my_fields->SiOI_density[field_idx1d]    = dsp[SpLUT::SiOI];
          my_fields->SiO2I_density[field_idx1d]   = dsp[SpLUT::SiO2I];
          my_fields->CH_density[field_idx1d]      = dsp[SpLUT::CH];
          my_fields->CH2_density[field_idx1d]     = dsp[SpLUT::CH2];
          my_fields->COII_density[field_idx1d]    = dsp[SpLUT::COII];
          my_fields->OII_density[field_idx1d]     = dsp[SpLUT::OII];
          my_fields->OHII_density[field_idx1d]    = dsp[SpLUT::OHII];
          my_fields->H2OII_density[field_idx1d]   = dsp[SpLUT::H2OII];
          my_fields->H3OII_density[field_idx1d]   = dsp[SpLUT::H3OII];
          my_fields->O2II_density[field_idx1d]    = dsp[SpLUT::O2II];
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
            if (my_chemistry->dust_species > 0)  {
              my_fields->Mg_density[field_idx1d]      = dsp[SpLUT::Mg];
            }
            if (my_chemistry->dust_species > 1)  {
              my_fields->Al_density[field_idx1d]      = dsp[SpLUT::Al];
              my_fields->S_density[field_idx1d]       = dsp[SpLUT::S];
              my_fields->Fe_density[field_idx1d]      = dsp[SpLUT::Fe];
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
          if (my_chemistry->dust_species > 0)  {
            my_fields->MgSiO3_dust_density[field_idx1d]  = dsp[SpLUT::MgSiO3_dust];
            my_fields->AC_dust_density[field_idx1d]      = dsp[SpLUT::AC_dust];
          }
          if (my_chemistry->dust_species > 1)  {
            my_fields->SiM_dust_density[field_idx1d]     = dsp[SpLUT::SiM_dust];
            my_fields->FeM_dust_density[field_idx1d]     = dsp[SpLUT::FeM_dust];
            my_fields->Mg2SiO4_dust_density[field_idx1d] = dsp[SpLUT::Mg2SiO4_dust];
            my_fields->Fe3O4_dust_density[field_idx1d]   = dsp[SpLUT::Fe3O4_dust];
            my_fields->SiO2_dust_density[field_idx1d]   = dsp[SpLUT::SiO2_dust];
            my_fields->MgO_dust_density[field_idx1d]     = dsp[SpLUT::MgO_dust];
            my_fields->FeS_dust_density[field_idx1d]     = dsp[SpLUT::FeS_dust];
            my_fields->Al2O3_dust_density[field_idx1d]   = dsp[SpLUT::Al2O3_dust];
          }
          if (my_chemistry->dust_species > 2)  {
            my_fields->ref_org_dust_density[field_idx1d]  = dsp[SpLUT::ref_org_dust];
            my_fields->vol_org_dust_density[field_idx1d]  = dsp[SpLUT::vol_org_dust];
            my_fields->H2O_ice_dust_density[field_idx1d]  = dsp[SpLUT::H2O_ice_dust];
          }
        }
      }

      e(i,j,k)     = dsp[i_eng-1];

      idsp.clear();
      vec.clear();
      mtrx = grackle::impl::View<double**>();
      mtrx_data_.clear();

    }
  }

  t_deriv::drop_MainScratchBuf(&main_scratch_buf);
  grackle::impl::drop_SpeciesCollection(&rhosp_dot);

}


} // namespace grackle::impl


#endif /* STEP_RATE_NEWTON_RAPHSON_HPP */
