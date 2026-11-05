//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements logic relevant to @ref step_rate_gauss_seidel
///
//===----------------------------------------------------------------------===//

#include "support/status_reporting.hpp"
#include "step_rate_gauss_seidel.hpp"

namespace GRIMPL_NAMESPACE_DECL {
namespace step_gs_detail {

DensUpdater::DensUpdater(const PartMap<SpKind>& sp_kind_map) {
  // perform a compile-time sanity-check to verify that the primordial species
  // participating in the network defined by primordial_chemistry==1 always
  // have the lowest 6 ids among all primordial species
  static_assert(((0 <= PrimordialSpLUT::e && PrimordialSpLUT::e < 6) &&
                 (0 <= PrimordialSpLUT::HI && PrimordialSpLUT::HI < 6) &&
                 (0 <= PrimordialSpLUT::HII && PrimordialSpLUT::HII < 6) &&
                 (0 <= PrimordialSpLUT::HeI && PrimordialSpLUT::HeI < 6) &&
                 (0 <= PrimordialSpLUT::HeII && PrimordialSpLUT::HeII < 6) &&
                 (0 <= PrimordialSpLUT::HeIII && PrimordialSpLUT::HeIII < 6)),
                "the ordering or primordial species ids changed unexpectedly");

  // infer the ids of species whose densities are updated before the electron
  // density (i.e. primordial species enabled by primordial_chemistry=1)
  IndexInterval1D prim_id_bounds = sp_kind_map.part_bounds(SpKind::PRIMORDIAL);
  int n_prim_sp = prim_id_bounds.stop - prim_id_bounds.start;
  GR_INTERNAL_REQUIRE(n_prim_sp >= 6, "sanity check failed!");
  species_id_e_ = prim_id_bounds.start + PrimordialSpLUT::e;
  species_id_HI_ = prim_id_bounds.start + PrimordialSpLUT::HI;
  species_id_HeIII_ = prim_id_bounds.start + PrimordialSpLUT::HeIII;

  spidarr_pre_elec_updates_[0] = species_id_HI_;
  spidarr_pre_elec_updates_[1] = prim_id_bounds.start + PrimordialSpLUT::HII;
  spidarr_pre_elec_updates_[2] = prim_id_bounds.start + PrimordialSpLUT::HeI;
  spidarr_pre_elec_updates_[3] = prim_id_bounds.start + PrimordialSpLUT::HeII;
  spidarr_pre_elec_updates_[4] = species_id_HeIII_;

  // record the information about the species whose densities are updated after
  // the electrron density
  IndexInterval1D remaining_prim_id_bounds{.start = prim_id_bounds.start + 6,
                                           .stop = prim_id_bounds.stop};
  remaining_id_info_pairs_[0] = {remaining_prim_id_bounds, false};

  // reminder: PartMap::part_bounds is designed such that when it's called
  //     for a kind that is not tracked by the kind_map, it will return an
  //     interval where start & stop are equal to each other (thus if it's
  //     used for the bounds of a for-loop, the loop will be skipped)
  remaining_id_info_pairs_[1] = {sp_kind_map.part_bounds(SpKind::METAL), true};
  remaining_id_info_pairs_[2] = {sp_kind_map.part_bounds(SpKind::DUST), true};
}

void DensUpdater::operator()(const double* dtit, IndexRange idx_range,
                             double* dedot_prev, double* HIdot_prev,
                             const gr_mask_type* itmask,
                             const gr_mask_type* itmask_metal,
                             const chemistry_data* my_chemistry,
                             SpeciesMultiView<gr_float> spdens,
                             View<const gr_float**> tmp_spdens) const {
  // Part 1: update buffers for each primordial species other than electrons
  //         that are enabled with primordial_chemistry >= 1
  {  // handle update to HIdot_prev
    gr_float* HI = spdens.contig1d_ptr(species_id_HI_);
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE) {
        double new_HI_val = tmp_spdens(species_id_HI_, i);
        // todo: break the gold standard and fix the next line
        HIdot_prev[i] =
            std::fabs(HI[i] - new_HI_val) / std::fmax((double)(dtit[i]), tiny8);

        if (std::isnan(new_HI_val)) {
          OMP_PRAGMA_CRITICAL {
            std::printf("HUGE HI! ::  %d %d %d %g\n", i, idx_range.j,
                        idx_range.k, new_HI_val);
          }
        }
      }
    }
  }

  for (int species_id : spidarr_pre_elec_updates_) {
    gr_float* ptr = spdens.contig1d_ptr(species_id);
    gr_float floor_val = (species_id == species_id_HeIII_)
                             ? static_cast<gr_float>(1e-5) * tiny_fortran_val
                             : tiny_fortran_val;
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if ((itmask[i] != MASK_FALSE)) {
        ptr[i] = std::fmax((gr_float)(tmp_spdens(species_id, i)), floor_val);
      }
    }
  }

  // Part 2: update the buffer for free electrons
  //
  // TODO: give some consideration to the case where all species have zero
  // density except for
  //   - species with neutral charges
  //   - Hminus & Dminus (H⁻ & D⁻).
  // Under the currrent charge conservation logic, the resulting electron
  // density would be negative in this scenario
  gr_float* de = spdens.contig1d_ptr(species_id_e_);
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      // temporarily store the electron density from the start of the current
      // subcycle.
      dedot_prev[i] = de[i];
    }
  }

  // actually override the electron density
  update_electron_densities(idx_range, itmask, itmask_metal, my_chemistry,
                            spdens);

  // update dedot_prev to store the time derivative of the electron density
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      dedot_prev[i] =
          std::fabs(de[i] - dedot_prev[i]) / std::fmax(dtit[i], tiny8);
    }
  }

  // Part 3: update buffers for all other species
  for (const spbounds_metalrelated_pair_& pair : remaining_id_info_pairs_) {
    const auto& [bounds, metal_related] = pair;

    for (int sp_id = bounds.start; sp_id < bounds.stop; sp_id++) {
      gr_float* ptr = spdens.contig1d_ptr(sp_id);
      if (metal_related) {
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          // at the moment, there are some cases where itmask holds
          // false and itmask_metal holds true
          if ((itmask[i] != MASK_FALSE) && (itmask_metal[i] != MASK_FALSE)) {
            ptr[i] =
                std::fmax((gr_float)(tmp_spdens(sp_id, i)), tiny_fortran_val);
          }
        }
      } else {
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (itmask[i] != MASK_FALSE) {
            ptr[i] =
                std::fmax((gr_float)(tmp_spdens(sp_id, i)), tiny_fortran_val);
          }
        }
      }
    }
  }
}

}  // namespace step_gs_detail
}  // namespace GRIMPL_NAMESPACE_DECL