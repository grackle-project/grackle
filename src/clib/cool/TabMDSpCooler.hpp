//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Declares @ref GRIMPL_NS::TabMDSpCooler
///
//===----------------------------------------------------------------------===//

#ifndef COOL_TABMDSPCOOLER_HPP
#define COOL_TABMDSPCOOLER_HPP

#include "../field_adaptor.hpp"
#include "../fortran_func_decls.h"  // gr_mask_type (todo: address this!)
#include "../internal_units.hpp"
#include "../interp_grid.hpp"
#include "../interpolate.hpp"
#include "../support/config.hpp"
#include "../support/index_helper.hpp"
#include "../support/View.hpp"

#include <vector>

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Encodes logic (and data) for computing species-specific cooling
///     contributions from multidimensional tables
struct TabMDSpCooler {
  // describes a single table
  struct SpCoolTablePack_ {
    double sp_mass_factor;
    int sp_id;
    bool neutral_H_ndens_ax_includes_molecular;
    const InterpGrid& interp_grid;
  };

  std::vector<SpCoolTablePack_> packs_;

  int pack_id_H2I_ = -1;
  int pack_id_HDI_ = -1;

private:
  // helper methods

  // does the heavy lifting!
  inline static void add_contrib_(
      double* edot, FortranView<const gr_float***> d, const gr_float* tgas,
      const gr_mask_type* itmask, const chemistry_data* my_chemistry,
      SpeciesMultiView<const double> sp_densities, const IndexRange& idx_range,
      const InternalGrUnits& internalu, double* scratch_buf,
      const SpCoolTablePack_* packs, int n_packs);

public:  // interface
  /// compute the number of scratch buffer elements given the number of
  /// elements considered at a given time
  static constexpr int num_scratch_buf_elements(int nelem_per_buf) noexcept {
    return 4 * nelem_per_buf;
  }

  /// @brief number of tracked rates for primordial species
  int n_prim_sp() const {
    return (static_cast<int>(pack_id_H2I_ >= 0) +
            static_cast<int>(pack_id_HDI_ >= 0));
  }

  /// @brief number of tracked rates for metal species
  int n_metal_sp() const { return packs_.size() - n_prim_sp(); }

  /// @brief Add cooling/heating contributions for H2I (if any)
  void add_H2I_contrib(double* edot, FortranView<const gr_float***> d,
                       const gr_float* tgas, const gr_mask_type* itmask,
                       const chemistry_data* my_chemistry,
                       SpeciesMultiView<const double> sp_densities,
                       const IndexRange& idx_range,
                       const InternalGrUnits& internalu,
                       double* scratch_buf) const {
    if (pack_id_H2I_ < 0) {
      return;
    }
    const SpCoolTablePack_& pack = packs_[pack_id_H2I_];
    add_contrib_(edot, d, tgas, itmask, my_chemistry, sp_densities, idx_range,
                 internalu, scratch_buf, &pack, 1);
  }

  /// @brief Add cooling/heating contributions for HDI (if any)
  void add_HDI_contrib(double* edot, FortranView<const gr_float***> d,
                       const gr_float* tgas, const gr_mask_type* itmask,
                       const chemistry_data* my_chemistry,
                       SpeciesMultiView<const double> sp_densities,
                       const IndexRange& idx_range,
                       const InternalGrUnits& internalu,
                       double* scratch_buf) const {
    if (pack_id_HDI_ < 0) {
      return;
    }
    const SpCoolTablePack_& pack = packs_[pack_id_HDI_];
    add_contrib_(edot, d, tgas, itmask, my_chemistry, sp_densities, idx_range,
                 internalu, scratch_buf, &pack, 1);
  }

  /// @brief Calculate cooling/heating for metal species
  void add_metal_sp_contrib(double* edot, FortranView<const gr_float***> d,
                            const gr_float* tgas, const gr_mask_type* itmask,
                            const chemistry_data* my_chemistry,
                            SpeciesMultiView<const double> sp_densities,
                            const IndexRange& idx_range,
                            const InternalGrUnits& internalu,
                            double* scratch_buf) const {
    const SpCoolTablePack_* packs = packs_.data() + n_prim_sp();
    add_contrib_(edot, d, tgas, itmask, my_chemistry, sp_densities, idx_range,
                 internalu, scratch_buf, packs, n_metal_sp());
  }
};

inline double interp_from_3D_grid(double input1, double input2, double input3,
                                  const GRIMPL_NS::InterpGrid& interp_grid) {
  return interpolate_3d(
      input1, input2, input3, interp_grid.props.dimension,
      interp_grid.props.parameters[0], interp_grid.props.parameter_spacing[0],
      interp_grid.props.parameters[1], interp_grid.props.parameter_spacing[1],
      interp_grid.props.parameters[2], interp_grid.props.parameter_spacing[2],
      interp_grid.props.data_size, interp_grid.data);
}

inline void TabMDSpCooler::add_contrib_(
    double* edot, FortranView<const gr_float***> d, const gr_float* tgas,
    const gr_mask_type* itmask, const chemistry_data* my_chemistry,
    SpeciesMultiView<const double> sp_densities, const IndexRange& idx_range,
    const InternalGrUnits& internalu, double* scratch_buf,
    const SpCoolTablePack_* packs, int n_packs) {
  if (n_packs == 0) {
    return;
  }

  int buf_len = sp_densities.extent(1);
  double* log10_T = scratch_buf;
  double* log10_H = scratch_buf + buf_len;
  double* log10_H2 = scratch_buf + 2 * buf_len;
  double* log10_dvdr = scratch_buf + 3 * buf_len;

  // calculate props reused for different interpolations
  double dom = internalu_calc_dom_(internalu);
  const double log10_Tcmb = (my_chemistry->cmb_temperature_floor == 1)
                                ? std::log10(internalu_calc_Tcmb_(internalu))
                                : 0.0;
  {
    const gr_float* HI = sp_densities.contig1d_ptr(SpLUT::HI);
    const gr_float* H2I = (my_chemistry->primordial_chemistry > 1)
                              ? sp_densities.contig1d_ptr(SpLUT::H2I)
                              : nullptr;

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        // TODO(gold-standard): reuse the pre-computed lnT values
        // -> we should either adjust the tables so that they interpolated with
        //    respect to lnT, or we should compute log10(T) = ln(T) / ln(10)
        log10_T[i] = std::log10(tgas[i]);
        log10_H[i] = std::log10(HI[i] * dom);
        if (my_chemistry->primordial_chemistry < 2) {
          log10_H2[i] = std::log10(HI[i] * dom);
        } else {
          log10_H2[i] = std::log10((HI[i] + H2I[i] / 2.0) * dom);
        }

        // compute log10(mass density in proper units of g/cm^3)
        double log10_rho =
            std::log10(d(i, idx_range.j, idx_range.k) * dom * constants::mH);

        // compute log10(dv/dr in units of km/s/cm)
        // -> From Chiaki & Wise (2019), approximate: dv/dr ≈ 1/(3 * t_ff)
        //    where t_ff = sqrt(3π/ (32 * G * ρ)) and ρ is mass density
        // -> let's break this down (since I originially found this confusing)
        //
        //    we start with
        //      (dv/dr) ≈ sqrt((32 * G * ρ) / (27π))
        //
        //    make the units explicit
        //      (dv/dr) / (1/s) = sqrt((32 * G_cgs * ρ / (g/cm³)) / (27π))
        //
        //    adjust so we can get our desired output units
        //      (dv/dr) / (km/s/cm) = 10⁻⁵ * sqrt((32 * G_cgs) / (27π))
        //                            * sqrt(ρ / (g/cm³))
        //
        //    to be clear (this was the source of my confusion) -8.799, is the
        //        value of:  `log10(10⁻⁵ * sqrt((32 * G_cgs) / (27π)))`
        //
        // TODO(gold-standard): explicitly use ``constants::GravConst``
        // -> for efficiency purposes, we want may to continue to precompute
        //    the term involving the gravitational constants
        log10_dvdr[i] = -8.79947961814 + 0.5 * log10_rho;
      }
    }
  }

  // perform each interpolation
  for (int pack_id = 0; pack_id < n_packs; pack_id++) {
    const SpCoolTablePack_& pack = packs[pack_id];
    const gr_float* cur_sp_dens = sp_densities.contig1d_ptr(pack.sp_id);
    const double* last_interp_ax_vals =
        (pack.neutral_H_ndens_ax_includes_molecular) ? log10_H2 : log10_H;

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        // nhat is `ndens_species / (dv/dr)` and has units of `s/(cm² * km)`
        double log10_ndens_sp =
            std::log10(cur_sp_dens[i] * dom / pack.sp_mass_factor);
        double log10_nhat = log10_ndens_sp - log10_dvdr[i];

        double log10_Linv = interp_from_3D_grid(
            log10_nhat, log10_T[i], last_interp_ax_vals[i], pack.interp_grid);
        double L = std::pow(10.0, -log10_Linv);

        double G = tiny8;
        if (my_chemistry->cmb_temperature_floor == 1) {
          double log10_Ginv = interp_from_3D_grid(
              log10_nhat, log10_Tcmb, last_interp_ax_vals[i], pack.interp_grid);
          G = std::pow(10.0, -log10_Ginv);
        }

        double edot_contrib =
            (G - L) / dom * cur_sp_dens[i] / pack.sp_mass_factor;
        if (std::isnan(edot_contrib)) {
          edot_contrib = 0.0;
        }
        edot[i] = edot[i] + edot_contrib;
      }
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // COOL_TABMDSPCOOLER_HPP