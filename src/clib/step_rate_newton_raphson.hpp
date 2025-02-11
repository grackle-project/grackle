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
#include "fortran_func_wrappers.hpp" // grackle::impl::fortran_wrapper::gaussj_g
#include "index_helper.h"
#include "internal_types.hpp"
#include "internal_units.h"
#include "utils-cpp.hpp"

#include "utils-field.hpp"
#include "time_deriv_0d.hpp"

namespace grackle::impl {

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

  // Density, energy and velocity fields fields

  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HI(my_fields->HI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(my_fields->HII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeI(my_fields->HeI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeII(my_fields->HeII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeIII(my_fields->HeIII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(my_fields->HM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(my_fields->H2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(my_fields->H2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DI(my_fields->DI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> DII(my_fields->DII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(my_fields->HDI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> e(my_fields->internal_energy, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  grackle::impl::View<gr_float***> DM(my_fields->DM_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(my_fields->HDII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(my_fields->HeHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CI(my_fields->CI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CII(my_fields->CII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO(my_fields->CO_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CO2(my_fields->CO2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OI(my_fields->OI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OH(my_fields->OH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2O(my_fields->H2O_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2(my_fields->O2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiI(my_fields->SiI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiOI(my_fields->SiOI_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> SiO2I(my_fields->SiO2I_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH(my_fields->CH_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> CH2(my_fields->CH2_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> COII(my_fields->COII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OII(my_fields->OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> OHII(my_fields->OHII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2OII(my_fields->H2OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H3OII(my_fields->H3OII_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> O2II(my_fields->O2II_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Mg(my_fields->Mg_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Al(my_fields->Al_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> S(my_fields->S_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> Fe(my_fields->Fe_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
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
  const int i_eng = 52;
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
      if (my_chemistry->primordial_chemistry > 0) { nsp = 6; }
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

      if ( my_chemistry->primordial_chemistry > 0 )  {
        dsp[ 1-1] = de(i,j,k);
        dsp[ 2-1] = HI(i,j,k);
        dsp[ 3-1] = HII(i,j,k);
        dsp[ 4-1] = HeI(i,j,k);
        dsp[ 5-1] = HeII(i,j,k);
        dsp[ 6-1] = HeIII(i,j,k);
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        dsp[ 7-1] = HM(i,j,k);
        dsp[ 8-1] = H2I(i,j,k);
        dsp[ 9-1] = H2II(i,j,k);
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        dsp[10-1] = DI(i,j,k);
        dsp[11-1] = DII(i,j,k);
        dsp[12-1] = HDI(i,j,k);
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        dsp[13-1] = DM(i,j,k);
        dsp[14-1] = HDII(i,j,k);
        dsp[15-1] = HeHII(i,j,k);
      }
      if ( itmask_metal[i] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          dsp[16-1] = CI(i,j,k);
          dsp[17-1] = CII(i,j,k);
          dsp[18-1] = CO(i,j,k);
          dsp[19-1] = CO2(i,j,k);
          dsp[20-1] = OI(i,j,k);
          dsp[21-1] = OH(i,j,k);
          dsp[22-1] = H2O(i,j,k);
          dsp[23-1] = O2(i,j,k);
          dsp[24-1] = SiI(i,j,k);
          dsp[25-1] = SiOI(i,j,k);
          dsp[26-1] = SiO2I(i,j,k);
          dsp[27-1] = CH(i,j,k);
          dsp[28-1] = CH2(i,j,k);
          dsp[29-1] = COII(i,j,k);
          dsp[30-1] = OII(i,j,k);
          dsp[31-1] = OHII(i,j,k);
          dsp[32-1] = H2OII(i,j,k);
          dsp[33-1] = H3OII(i,j,k);
          dsp[34-1] = O2II(i,j,k);
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0)  {
              dsp[35-1] = Mg(i,j,k);
            }
            if (my_chemistry->dust_species > 1)  {
              dsp[36-1] = Al(i,j,k);
              dsp[37-1] = S(i,j,k);
              dsp[38-1] = Fe(i,j,k);
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            dsp[39-1] = MgSiO3(i,j,k);
            dsp[40-1] = AC(i,j,k);
          }
          if (my_chemistry->dust_species > 1)  {
            dsp[41-1] = SiM(i,j,k);
            dsp[42-1] = FeM(i,j,k);
            dsp[43-1] = Mg2SiO4(i,j,k);
            dsp[44-1] = Fe3O4(i,j,k);
            dsp[45-1] = SiO2D(i,j,k);
            dsp[46-1] = MgO(i,j,k);
            dsp[47-1] = FeS(i,j,k);
            dsp[48-1] = Al2O3(i,j,k);
          }
          if (my_chemistry->dust_species > 2)  {
            dsp[49-1] = reforg(i,j,k);
            dsp[50-1] = volorg(i,j,k);
            dsp[51-1] = H2Oice(i,j,k);
          }
        }
      }
      dsp[i_eng-1] = e(i,j,k);

      id = 0;
      if (my_chemistry->primordial_chemistry > 0)  {
        for (isp = 1; isp<=(6); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 1)  {
        for (isp = 7; isp<=(9); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 2)  {
        for (isp = 10; isp<=(12); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (my_chemistry->primordial_chemistry > 3)  {
        for (isp = 13; isp<=(15); isp++) {
          id = id + 1;
          idsp[id-1] = isp;
        }
      }
      if (itmask_metal[i] != MASK_FALSE)  {
        if (my_chemistry->metal_chemistry == 1)  {
          for (isp = 16; isp<=(34); isp++) {
            id = id + 1;
            idsp[id-1] = isp;
          }
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
            if (my_chemistry->dust_species > 0)  {
              for (isp = 35; isp<=(35); isp++) {
                id = id + 1;
                idsp[id-1] = isp;
              }
            }
            if (my_chemistry->dust_species > 1)  {
              for (isp = 36; isp<=(38); isp++) {
                id = id + 1;
                idsp[id-1] = isp;
              }
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1) )  {
          if (my_chemistry->dust_species > 0)  {
            for (isp = 39; isp<=(40); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
          if (my_chemistry->dust_species > 1)  {
            for (isp = 41; isp<=(48); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
          if (my_chemistry->dust_species > 2)  {
            for (isp = 49; isp<=(51); isp++) {
              id = id + 1;
              idsp[id-1] = isp;
            }
          }
        }
      }
      if ( imp_eng[i] ==1 )  {
        id = id + 1;
        idsp[id-1] = i_eng;
      }

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

      // before we merge the gen2024 branch into the main grackle branch we
      // need to resolve the following issue (or at least exit gracefully)
      // -> the easiest short-term solution is to temporarily overwrite the
      //    value of iter
      // -> Alternatively, we could try to reuse the existing value of tgasold.
      //    Doing that requires some care; we need to use meaningful values
      //    of tgasold when we iteratively call lookup_cool_rates0d for the
      //    sake of estimating elements in the jacobian matrix.
      GRIMPL_REQUIRE(
        ((iter == 1) || (imp_eng[i] != 1)),
        "there is a logical issue in the original lookup_cool_rates0d "
        "Fortran routine when (iter != 1) AND (imp_eng == 1).\n"
        " -> as it was originally written the routine, allocates storage\n"
        "    for the tgasold buffer (on the stack) and leaves the value\n"
        "    unitialized.\n"
        " -> this is a problem when (iter != 1) and (imp_eng == 1) because\n"
        "    cool1d_multi_g will try to make use of the unitialized variable\n"
        " -> this isn't a problem when (iter == 1) and (imp_eng == 1) since\n"
        "    cool1d_multi_g knows it needs to initialize tgasold. It isn't\n"
        "    a problem when (imp_eng != 1) since cool1d_multi_g isn't called"
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
          t_deriv::derivatives(dt_FIXME, dsp.data(), dspdot.data(), pack);

          for (jsp = 1; jsp<=(nsp); jsp++) {
            dspj = eps * dsp[idsp[jsp-1]-1];
            for (isp = 1; isp<=(nsp); isp++) {
              if(isp == jsp)  {
                dsp1[idsp[isp-1]-1] = dsp[idsp[isp-1]-1] + dspj;
              } else {
                dsp1[idsp[isp-1]-1] = dsp[idsp[isp-1]-1];
              }
            }

            t_deriv::derivatives(dt_FIXME, dsp1.data(), dspdot1.data(), pack);

            for (isp = 1; isp<=(nsp); isp++) {
              if ( (dsp[idsp[isp-1]-1]==0.e0)
               &&  (dspdot1[idsp[isp-1]-1]
               ==  dspdot[idsp[isp-1]-1]) )  {
                der(idsp[isp-1]-1,idsp[jsp-1]-1) = 0.e0;
              } else {
                der(idsp[isp-1]-1,idsp[jsp-1]-1) =
                   (dspdot1[idsp[isp-1]-1]
                   - dspdot[idsp[isp-1]-1]) / dspj;
              }
            }

          }

          for (isp = 1; isp<=(nsp); isp++) {
            for (jsp = 1; jsp<=(nsp); jsp++) {
              if(isp == jsp)  {
                mtrx(isp-1,jsp-1) = 1.e0 - dtit[i]
                   * der(idsp[isp-1]-1,idsp[jsp-1]-1);
              } else {
                mtrx(isp-1,jsp-1) =      - dtit[i]
                   * der(idsp[isp-1]-1,idsp[jsp-1]-1);
              }
            }
          }

          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = dspdot[idsp[isp-1]-1] * dtit[i]
                   - ddsp[idsp[isp-1]-1];
          }

          // to get more accuracy
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]/d(i,j,k);
          }

          ierror = f_wrap::gaussj_g(nsp, mtrx.data(), vec.data());
          if(ierror == 1)  {
            goto label_9998;
          }

          // multiply with density again
          for (isp = 1; isp<=(nsp); isp++) {
            vec[isp-1] = vec[isp-1]*d(i,j,k);
          }

          for (isp = 1; isp<=(nsp); isp++) {
            ddsp[idsp[isp-1]-1] = ddsp[idsp[isp-1]-1] + vec[isp-1];
            dsp[idsp[isp-1]-1]  = dsp[idsp[isp-1]-1]  + vec[isp-1];
          }

          if (imp_eng[i] == 1)  {
            if( (my_chemistry->primordial_chemistry > 0)  &&  (my_chemistry->with_radiative_cooling == 1) )  {
              for (isp = 1; isp<=(nsp); isp++) {
                if ( (dsp[idsp[isp-1]-1] != dsp[idsp[isp-1]-1])
                 ||  (dsp[idsp[isp-1]-1] <= 0.) )  {
                  ierror = 1;
                  goto label_9997;
                }
              }
            }
          }

          err_max = 0.e0;
          for (isp = 1; isp<=(nsp); isp++) {
            if(dsp[idsp[isp-1]-1] > tiny8)  {
              err = grackle::impl::dabs(vec[isp-1] / dsp[idsp[isp-1]-1]);
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
            if ( (dsp[idsp[isp-1]-1] != dsp[idsp[isp-1]-1])
             ||  (dsp[idsp[isp-1]-1] <= 0.) )  {
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

      if ( my_chemistry->primordial_chemistry > 0 )  {
        de(i,j,k)      = dsp[ 1-1];
        HI(i,j,k)      = dsp[ 2-1];
        HII(i,j,k)     = dsp[ 3-1];
        HeI(i,j,k)     = dsp[ 4-1];
        HeII(i,j,k)    = dsp[ 5-1];
        HeIII(i,j,k)   = dsp[ 6-1];
      }
      if ( my_chemistry->primordial_chemistry > 1 )  {
        HM(i,j,k)      = dsp[ 7-1];
        H2I(i,j,k)     = dsp[ 8-1];
        H2II(i,j,k)    = dsp[ 9-1];
      }
      if ( my_chemistry->primordial_chemistry > 2 )  {
        DI(i,j,k)      = dsp[10-1];
        DII(i,j,k)     = dsp[11-1];
        HDI(i,j,k)     = dsp[12-1];
      }
      if ( my_chemistry->primordial_chemistry > 3 )  {
        DM(i,j,k)      = dsp[13-1];
        HDII(i,j,k)    = dsp[14-1];
        HeHII(i,j,k)   = dsp[15-1];
      }
      if ( itmask_metal[i] != MASK_FALSE )  {
        if ( my_chemistry->metal_chemistry == 1 )  {
          CI(i,j,k)      = dsp[16-1];
          CII(i,j,k)     = dsp[17-1];
          CO(i,j,k)      = dsp[18-1];
          CO2(i,j,k)     = dsp[19-1];
          OI(i,j,k)      = dsp[20-1];
          OH(i,j,k)      = dsp[21-1];
          H2O(i,j,k)     = dsp[22-1];
          O2(i,j,k)      = dsp[23-1];
          SiI(i,j,k)     = dsp[24-1];
          SiOI(i,j,k)    = dsp[25-1];
          SiO2I(i,j,k)   = dsp[26-1];
          CH(i,j,k)      = dsp[27-1];
          CH2(i,j,k)     = dsp[28-1];
          COII(i,j,k)    = dsp[29-1];
          OII(i,j,k)     = dsp[30-1];
          OHII(i,j,k)    = dsp[31-1];
          H2OII(i,j,k)   = dsp[32-1];
          H3OII(i,j,k)   = dsp[33-1];
          O2II(i,j,k)    = dsp[34-1];
          if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
            if (my_chemistry->dust_species > 0)  {
              Mg(i,j,k)      = dsp[35-1];
            }
            if (my_chemistry->dust_species > 1)  {
              Al(i,j,k)      = dsp[36-1];
              S(i,j,k)       = dsp[37-1];
              Fe(i,j,k)      = dsp[38-1];
            }
          }
        }
        if ( ( my_chemistry->grain_growth == 1 )  ||  ( my_chemistry->dust_sublimation == 1 ) )  {
          if (my_chemistry->dust_species > 0)  {
            MgSiO3(i,j,k)  = dsp[39-1];
            AC(i,j,k)      = dsp[40-1];
          }
          if (my_chemistry->dust_species > 1)  {
            SiM(i,j,k)     = dsp[41-1];
            FeM(i,j,k)     = dsp[42-1];
            Mg2SiO4(i,j,k) = dsp[43-1];
            Fe3O4(i,j,k)   = dsp[44-1];
            SiO2D(i,j,k)   = dsp[45-1];
            MgO(i,j,k)     = dsp[46-1];
            FeS(i,j,k)     = dsp[47-1];
            Al2O3(i,j,k)   = dsp[48-1];
          }
          if (my_chemistry->dust_species > 2)  {
            reforg(i,j,k)  = dsp[49-1];
            volorg(i,j,k)  = dsp[50-1];
            H2Oice(i,j,k)  = dsp[51-1];
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

}


} // namespace grackle::impl


#endif /* STEP_RATE_NEWTON_RAPHSON_HPP */
