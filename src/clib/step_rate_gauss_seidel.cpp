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

}  // namespace step_gs_detail
}  // namespace GRIMPL_NAMESPACE_DECL