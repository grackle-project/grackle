//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the calc_grain_size_increment_species_1d function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// calc_grain_size_increment_species_1d function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "calc_grain_size_increment_species_1d.hpp"

void grackle::impl::calc_grain_size_increment_species_1d(
    int igrgr, const gr_mask_type* itmask, int SN0_N, int in, int jn, int kn,
    IndexRange idx_range, gr_float* density_data, int nSN,
    const gr_float* dsp_data_, gr_float* SN_metal_data, double* SN_fsp,
    double* SN_r0sp_data, double ssp, double* sgsp, double* alsp_data,
    int* gr_N, int gr_Size, double* SN_kp0sp_data) {
  // input
  int iSN;

  grackle::impl::View<gr_float***> d(density_data, in, jn, kn);
  grackle::impl::View<const gr_float***> dsp(dsp_data_, in, jn, kn);
  grackle::impl::View<gr_float**> SN_metal(SN_metal_data, in, SN0_N);

  // table
  grackle::impl::View<double**> SN_r0sp(SN_r0sp_data, 3, SN0_N);

  // opacity table
  grackle::impl::View<double**> SN_kp0sp(SN_kp0sp_data, gr_Size, SN0_N);

  // output
  grackle::impl::View<double**> alsp(alsp_data, gr_N[1], in);

  // local
  int i;
  double coef0, coef1, coef2, coef3;
  double dsp_inject_sum;
  double SN_sgsp, SN_kpsp;
  std::vector<double> SN_dsp0(SN0_N);
  std::vector<double> SN_nsp0(SN0_N);
  std::vector<double> drsp(in);
  const double pi_local_var = pi_fortran_val;
  int iTd, iTd0;

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      // Step 1: compute the total mass density of the current grain species
      //         that was injected (by summing the amounts injected by each
      //         injection pathway)
      for (iSN = 0; iSN < nSN; iSN++) {
        if (SN_fsp[iSN] > 0.e0) {
          SN_dsp0[iSN] = SN_fsp[iSN] * SN_metal(i, iSN);
        }
      }

      // Step 2: Compute the size increment for the grain species

      // Let's go into detail:
      // - todo: maybe we move most of this description into the docstring
      //   and possibly move some of it into the narrative documentation
      // - let's define some notation:
      //    - for some time t, let φ(r,t) denote the grain species's
      //      differential size distribution at time t, ρ(t) denote the grain
      //      species's mass density, and n(t) denote it number density.
      //    - we define <rᵖ>(t) as the `p`th order moment of `φ(r,t)`. In
      //      other words, <rᵖ>ⱼ = ∫ rᵖ φⱼ dr, where the integral is taken
      //      from 0 to infinity
      //    - for additional context:
      //      - φ(r,t) is normalized such that the p=0 moment is 1
      //      - n(t)*φ(r,t) specifies the number density of the grain with
      //        radii between r and r+dr
      //      - given the mass density of a single grain, ζ, the variables
      //        are related via ρ(t) = (4π/3) ζ n(t) <r³>(t).
      //    - we define Φⱼ(r) as the grain species's initial differential
      //      size distribution when distributed by injection pathway j.
      //    - for notational convenience, we define <rᵖ>ⱼ as the `p`th order
      //      moment of `Φⱼ(r)`
      //    - NOTE: we don't explicitly track φ(r,t) or Φⱼ(r), instead we
      //      track moments
      // - The underlying model assumes that all injection pathways inject
      //   the grain species at a single time t_inj.
      //   - In other words:  φ(r,t_inj) = (∑ⱼ nⱼ Φⱼ(r)) / (∑ⱼ nⱼ)
      //   - Consequently:
      //         <rᵖ>(t_inj) = (∑ⱼ nⱼ <rᵖ>ⱼ) / (∑ⱼ nⱼ)
      //         OR
      //         <rᵖ>(t_inj) = (∑ⱼ <rᵖ>ⱼ ρⱼ/<r³>ⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ)
      // - Furthermore we assume that:
      //   - grains can only be created via injection (i.e. the number
      //     density of the grain species is frozen)
      //   - φ(r,t) can't be deformed. However it can be translated (which
      //     is exactly what happens when grains undergo growth).
      // - Thus, when we model growth, φ(r,t) = φ(r - δr(t), t_inj) where
      //   we call δr(t) the "size increment"

      // Anybody reading this might notice that several limitations to this
      // model. This is discussed in detail in GH Issue 444 (the description
      // should be made part of the narrative docs

      if (igrgr == 0) {
        // the model is constructed such that δr(t) = 0 if we aren't
        // modelling grain growth

        drsp[i] = 0.e0;

      } else {
        // here we calculate the current grain size increment, δr(t), using
        // conservation of grain number (i.e. grains are only created via
        // injection pathways)
        // - recall that to model growth we effectively assume that each
        //   occurrence of a grain species was injected at a single time, t_inj
        // - conservation of grain number tells us that
        //   <r³>(t) = <r³>(t_inj) * ρ(t) / ρ(t_inj),
        //   where ρ(t_inj) = ∑ⱼ ρⱼ
        // - Equation A3 of Chiaki & Wise 2019 indicates that the grain-size
        //   increment is the root of the following cubic equation
        //    [<r³>(t_inj) - <r³>(t)] + [3 * <r²>(t_inj)] * δr(t)
        //       + [3 * <r¹>(t_inj)] * δr²(t) + δr³(t) = 0
        // - Let's rewrite the expression (<r³>(t_inj) - <r³>(t)). First,
        //   let's plug in the equation from conservation of grain number:
        //     <r³>(t_inj) - <r³>(t) = <r³>(t_inj) * [ 1 - (ρ(t)/∑ⱼρⱼ)]
        //   Now, let's plug in <r³>(t_inj) = (∑ⱼρⱼ) / (∑ⱼ ρⱼ/<r³>ⱼ), which
        //   comes from the generic formula for <rᵖ>(t_inj) defined earlier:
        //     <r³>(t_inj) - <r³>(t) = [(∑ⱼρⱼ) - ρ(t)] / (∑ⱼ ρⱼ/<r³>ⱼ)
        // - If we plug this into the preceding cubic equation and multiply the
        //   resulting equation by (∑ⱼ ρⱼ/<r³>ⱼ), we get:
        //       coef0 + coef1 * δr(t) + coef2 * δr²(t) + coef3 * δr³(t) = 0
        //   where:
        //       - coef0 = (∑ⱼ ρⱼ) - ρ(t)
        //       - coef1 = ∑ⱼ 3 <r²>ⱼ ρⱼ/<r³>ⱼ
        //       - coef2 = ∑ⱼ 3 <r¹>ⱼ ρⱼ/<r³>ⱼ
        //       - coef3 = ∑ⱼ ρⱼ/<r³>ⱼ

        // Let's compute these coefficients
        dsp_inject_sum = 0.e0;
        // We'll use dsp_inject_sum to compute coef0
        coef1 = 0.e0;
        coef2 = 0.e0;
        coef3 = 0.e0;

        // Loop over each injection pathway
        for (iSN = 0; iSN < nSN; iSN++) {
          if (SN_fsp[iSN] > 0.e0) {
            // Calculate 4πζnⱼ/3 = ρⱼ/<r³>ⱼ
            // -> recall: that ζ is the mass density of a single grain of the
            //    current grain species (i.e. it's a constant)
            // -> TODO: its very confusing that we store the result within
            // SN_nsp0
            //    since that variable is later reused to directly hold nⱼ
            SN_nsp0[iSN] = SN_dsp0[iSN] / SN_r0sp(2, iSN);

            dsp_inject_sum = dsp_inject_sum + SN_dsp0[iSN];
            coef1 = coef1 + 3.e0 * SN_nsp0[iSN] * SN_r0sp(1, iSN);
            coef2 = coef2 + 3.e0 * SN_nsp0[iSN] * SN_r0sp(0, iSN);
            coef3 = coef3 + SN_nsp0[iSN];
          }
        }

        coef0 = dsp_inject_sum - dsp(i, idx_range.j, idx_range.k);

        // Let's actually solve for the δr(t), the root of the cubic equation
        // - before we do that, we divide both sides of the equation by the
        //   value of coef3. This ensures that the rescaled value of coef3 is
        //   exactly 1 (a requirement of the solve_cubic_equation routine)
        // - todo: consider not modifying the variables in-place or performing
        //   this operation within the solve_cubic_equation
        coef0 = coef0 / coef3;
        coef1 = coef1 / coef3;
        coef2 = coef2 / coef3;

        FORTRAN_NAME(solve_cubic_equation)(&coef2, &coef1, &coef0, &drsp[i]);
        // TODO: to be removed after fixing the numerical issue with i=0 and
        // idx_range.j=idx_range.k=0
        if (i == 0 && idx_range.j == 0 && idx_range.k == 0) {
          drsp[i] = 0.e0;
        }

        drsp[i] = std::fmax(drsp[i], 0.e0);
      }

      // Step 3: calculate number density (code_density / g)

      for (iSN = 0; iSN < nSN; iSN++) {
        if (SN_fsp[iSN] > 0.e0) {
          SN_nsp0[iSN] = SN_dsp0[iSN] /
                         (4.e0 * pi_local_var / 3.e0 * ssp * SN_r0sp(2, iSN));
        } else {
          SN_nsp0[iSN] = 0.e0;
        }
      }

      // Step 4: calculate geometrical cross-section per unit gas mass
      // -> units of cm^2/g
      sgsp[i] = 0.e0;
      for (iSN = 0; iSN < nSN; iSN++) {
        if (SN_fsp[iSN] > 0.e0) {
          SN_sgsp = pi_local_var *
                    (SN_r0sp(1, iSN) + 2.e0 * SN_r0sp(0, iSN) * drsp[i] +
                     std::pow(drsp[i], 2));
        } else {
          SN_sgsp = 0.e0;
        }
        sgsp[i] = sgsp[i] + SN_nsp0[iSN] * SN_sgsp;
      }
      sgsp[i] = sgsp[i] / d(i, idx_range.j, idx_range.k);

      // Step 5: calculate optical opacity related quantities
      // -> we are effectively constructing a 1d table of values, at various
      //    possible grain temperature for each injection pathway (in other
      //    words, its a 2D table)
      // -> I'm pretty sure that the values we are the opacity (per unit gas
      //    mass)
      // -> I'm pretty confident that the units are independent of code units
      //    (I think the units are cm^2/g)
      for (iTd = 0; iTd < (gr_N[1]); iTd++) {
        iTd0 = iTd * gr_N[0];
        alsp(iTd, i) = 0.e0;
        for (iSN = 0; iSN < nSN; iSN++) {
          if (SN_fsp[iSN] > 0.e0) {
            SN_kpsp = 4.e0 * pi_local_var / 3.e0 * ssp *
                      (SN_kp0sp(iTd0 + 3, iSN) +
                       3.e0 * SN_kp0sp(iTd0 + 2, iSN) * drsp[i] +
                       3.e0 * SN_kp0sp(iTd0 + 1, iSN) * std::pow(drsp[i], 2) +
                       SN_kp0sp(iTd0 + 0, iSN) * std::pow(drsp[i], 3));
          } else {
            SN_kpsp = 0.e0;
          }
          alsp(iTd, i) = alsp(iTd, i) + SN_nsp0[iSN] * SN_kpsp;
        }
        alsp(iTd, i) = alsp(iTd, i) / d(i, idx_range.j, idx_range.k);
      }
    }
  }

  return;
}
