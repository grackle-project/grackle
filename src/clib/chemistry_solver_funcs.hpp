//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// @brief Defines chemistry reaction related functions invoked by the
///     grackle solver in order to integrate the species densities over time.
///
/// The premise is for this file to hold functions used during the calculation
/// of species density derivatives during integration and to decouple this
/// logic from the integration to the greatest extent possible. This includes:
/// - `lookup_cool_rates1d_g`, when we eventually transcribe it
/// - the logic that is duplicated between `step_rate_g` and
///   `species_density_derivatives_0d`.
///   - the ordering of calculations required by `step_rate_g` proabably can't
///     be decoupled from the derivative calculation for primoridial species
///   - it may also make sense to further divide logic by the kinds of species
///     that are affected (e.g. primordial vs grains)
///
//===----------------------------------------------------------------------===//

#ifndef CHEMISTRY_SOLVER_FUNCS_HPP
#define CHEMISTRY_SOLVER_FUNCS_HPP

#include <climits>

#include "chem_model/infer_species.hpp"
#include "grackle.h"
#include "field_adaptor.hpp"
#include "fortran_func_decls.h"  // gr_mask_type
#include "full_rxn_rate_buf.hpp"
#include "support/index_helper.hpp"
#include "LUT.hpp"
#include "support/status_reporting.hpp"
#include "support/View.hpp"

namespace grackle::impl::chemistry {

/// @brief Converts the relative species id to the global species id
///
/// This is a tool to help us replace `SpLUT`. Given an object `lut`, you
/// can replace occurrences of `SpLUT::<name>` with one of:
/// - `spid.prim(PrimordialSpLUT::<name>)`
/// - `spid.metal(MetalSpLUT::<name>)`
/// - `spid.dust(DustSpLUT::<name>)`
///
/// An important property of this type is that even when grackle isn't
/// configured to track a given species, each of the above method calls will
/// fallback to a valid species id.
/// - in a vaccuum, this property looks and acts a lot like bug (if grackle is
///   misconfigured, we could accidentally use the wrong species)
/// - however, this property is required to facillitate rewrites of the
///   `chemistry::species_density_*` without introducing undefined behavior or
///   requiring significant immediate structural changes
/// The primary reason this property is required is the historical choice to
/// eagerly load a species density buffer at the start of the function (even
/// if the species isn't actually used)
///
/// Longer Term Outlook
/// ===================
/// Longer term, the goal is to rewrite the `chemistry::species_density_*`
/// functions in order to totally remove this class
class SpIdAdaptor{
  // tracks properties about an interval of species ids that are all a given
  // kind of species (e.g. primordial, metal, or dust)
  struct IntervalProps{
    int id_offset;
    int id_count;
  };

  int fallback_id_;
  IntervalProps prim_;
  IntervalProps metal_;
  IntervalProps dust_;

  [[gnu::always_inline]] int calc_(const IntervalProps& kind_id_props,
                                   int kind_id) const noexcept {
    // let's skip the check that kind_id is not-negative
    bool valid = kind_id < kind_id_props.id_count;
    // we could make the following branchless (unclear if that matters)
    return (valid) ? kind_id_props.id_offset + kind_id : fallback_id_;
  };

public:
  // define the main methods:
  [[gnu::always_inline]] int prim(int id) const {return calc_(prim_, id); }
  [[gnu::always_inline]] int metal(int id) const {return calc_(metal_, id); }
  [[gnu::always_inline]] int dust(int id) const { return calc_(dust_, id); }

  SpIdAdaptor() = delete;

  /// @brief Primary Constructor
  explicit SpIdAdaptor(PartMap<SpKind> sp_kind_map) {
    int min_id = INT_MAX;

    SpKind kinds[3] = {SpKind::PRIMORDIAL, SpKind::METAL, SpKind::DUST};
    for (SpKind kind: kinds) {
      IndexInterval1D b = sp_kind_map.part_bounds(kind);
      IntervalProps props{.id_offset = b.start, .id_count = b.stop - b.start};
      if (kind == SpKind::PRIMORDIAL) {
        prim_ = props;
      } else if (kind == SpKind::METAL) {
        metal_ = props;
      } else {
        dust_ = props;
      }
      min_id = (b.start >= 0 && b.start < min_id) ? b.start : min_id;
    }
    fallback_id_ = min_id;

    if (min_id == INT_MAX) { // this should be impossible
      GR_INTERNAL_ERROR("sp_kind_map contains no species");
    }
  }
};

/// Perform an implicit Gauss-Seidel sweep of a backward-Euler time integrator
/// to advance the rate equations by one (sub-)cycle (for each index in the
/// index-range that is selected by the given itmask)
///
/// @param[out] out_spdens stores the updated species densities for the
///     @p idx_range.
/// @param[in] idx_range Specifies the current index-range
/// @param[in] dtit Specifies the timestep of the current sub-cycle for each
///     index in @p idx_range.
/// @param[in] anydust Indicates whether we are modelling dust
/// @param[in] rhoH Indicates the mass density of all Hydrogen
/// @param[in] itmask The general iteration mask for @p idx_range.
/// @param[in] itmask_metal The iteration mask @p idx_range that specifies
///     where we should account for metal-species chemistry and grain-species
///     chemistry.
/// @param[in] my_chemistry Provides various runtime parameters (we probably
///     don't need to pass the whole thing)
/// @param[in] sp_densities Specifies the densities of the various species
///     that Grackle evolves in a format that allows the values to be accessed
///     with the index lookup table.
/// @param[in] rxn_rate_buf specifies the precomputed rxn rates (depends on
///     local physical conditions)
/// @param[in]  sp_kind_map associates intervals of species ids with the kind
///     with the appropriate @ref SpKind
///
/// Refactoring Goals
/// -----------------
///
/// Right now the implementation of every rate is **very** manual, which makes
/// it very easy to forget to add a rate.
///
/// The shorter-term goal is to deduplicate as much as possible between this
/// function & @ref species_density_derivatives_0d.
/// Some important considerations:
/// - it will be easiest to deduplicate the metal-chemistry and grain-species
///   growth rates (we only need to introduce some light usage of templates)
/// - deduplicating other logic will involve additional templates
/// - introducing light-weight templates to deduplicate the metal-chemistry
///   and grain-species growth rate logic is well worth the cost to improve
///   the maintenance burden. Deduplicating the other logic is probably also
///   worth the cost (but it's more debatable)
///
/// The longer-term goal is to overhaul this function.
/// - In the longer term, it would make more sense to use more of a "table
///   based approach" in the regime where we have lots of chemistry.
/// - We elaborate a little more down below:
///   - we already have a table of 1d rate buffers: @p rxn_rate_buf.
///   - we might also track a table (with matching indices) where we track the
///     indexes associated with each reactant/product of a rate as well as
///     stoichiometric coefficients.
///   - Furthermore, we could envision tracking rate-lists for each species of
///     all creational reactions and destructive reactions.
///   - If we have all of this, then we could just iterate over each species'
///     reactions and use the tables to dynamically access reactant densities
///     by index. Essentially, we would dynamically build up the creational
///     & destructive coefficients (``acoef`` & ``scoef``) for a full row and
///     use the values to then perform the update.
///   - in practice things would be a little more complex since we have to deal
///     dissociation rates, but it's still **very** doable
///   - once this infrastructure is in place, we could replace
///     @ref chemistry::species_density_derivatives_0d with logic that iterates
///     through rates in a very different manner and directly compute partial
///     derivatives in a much more efficient manner (currently, it's used with
///     finite differences for computing the partial derivatives)
/// - Importantly, the table-based approach is *probably* slower than the
///   manual, hand-coded approach in the limit of a small chemical network.
///   I suspect that we'll preserve the hand-coded approach for
///   primordial_chemistry == 1, probably primordial_chemistry == 2, and
///   maybe even primordial_chemistry == 3 (but, we'll obviously need to test
///   performance).
inline void species_density_updates_gauss_seidel(
  View<double**> out_spdens, IndexRange idx_range,
  const double* dtit, gr_mask_type anydust, const double* rhoH,
  const gr_mask_type* itmask, const gr_mask_type* itmask_metal,
  const chemistry_data* my_chemistry,
  SpeciesMultiView<const gr_float> spdens,
  const FullRxnRateBuf rxn_rate_buf, const PartMap<SpKind>& sp_kind_map
)
{
  const SpIdAdaptor spid(sp_kind_map);


  // Get pointers to the various species fields
  // ------------------------------------------
  const gr_float* de = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::e));
  const gr_float* HI = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HI));
  const gr_float* HII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HII));
  const gr_float* HeI = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HeI));
  const gr_float* HeII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HeII));
  const gr_float* HeIII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HeIII));
  const gr_float* HM = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HM));
  const gr_float* H2I = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::H2I));
  const gr_float* H2II = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::H2II));
  const gr_float* DI = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::DI));
  const gr_float* DII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::DII));
  const gr_float* HDI = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HDI));
  const gr_float* DM = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::DM));
  const gr_float* HDII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HDII));
  const gr_float* HeHII = spdens.contig1d_ptr(spid.prim(PrimordialSpLUT::HeHII));
  const gr_float* CI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::CI));
  const gr_float* CII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::CII));
  const gr_float* COI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::COI));
  const gr_float* CO2I = spdens.contig1d_ptr(spid.metal(MetalSpLUT::CO2I));
  const gr_float* OI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::OI));
  const gr_float* OHI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::OHI));
  const gr_float* H2OI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::H2OI));
  const gr_float* O2I = spdens.contig1d_ptr(spid.metal(MetalSpLUT::O2I));
  const gr_float* SiI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::SiI));
  const gr_float* SiOI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::SiOI));
  const gr_float* SiO2I = spdens.contig1d_ptr(spid.metal(MetalSpLUT::SiO2I));
  const gr_float* CHI = spdens.contig1d_ptr(spid.metal(MetalSpLUT::CHI));
  const gr_float* CH2I = spdens.contig1d_ptr(spid.metal(MetalSpLUT::CH2I));
  const gr_float* COII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::COII));
  const gr_float* OII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::OII));
  const gr_float* OHII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::OHII));
  const gr_float* H2OII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::H2OII));
  const gr_float* H3OII = spdens.contig1d_ptr(spid.metal(MetalSpLUT::H3OII));
  const gr_float* O2II = spdens.contig1d_ptr(spid.metal(MetalSpLUT::O2II));
  const gr_float* Mg = spdens.contig1d_ptr(spid.metal(MetalSpLUT::Mg));
  const gr_float* Al = spdens.contig1d_ptr(spid.metal(MetalSpLUT::Al));
  const gr_float* S = spdens.contig1d_ptr(spid.metal(MetalSpLUT::S));
  const gr_float* Fe = spdens.contig1d_ptr(spid.metal(MetalSpLUT::Fe));
  const gr_float* SiM = spdens.contig1d_ptr(spid.dust(DustSpLUT::SiM_dust));
  const gr_float* FeM = spdens.contig1d_ptr(spid.dust(DustSpLUT::FeM_dust));
  const gr_float* Mg2SiO4 = spdens.contig1d_ptr(spid.dust(DustSpLUT::Mg2SiO4_dust));
  const gr_float* MgSiO3 = spdens.contig1d_ptr(spid.dust(DustSpLUT::MgSiO3_dust));
  const gr_float* Fe3O4 = spdens.contig1d_ptr(spid.dust(DustSpLUT::Fe3O4_dust));
  const gr_float* AC = spdens.contig1d_ptr(spid.dust(DustSpLUT::AC_dust));
  const gr_float* SiO2D = spdens.contig1d_ptr(spid.dust(DustSpLUT::SiO2_dust));
  const gr_float* MgO = spdens.contig1d_ptr(spid.dust(DustSpLUT::MgO_dust));
  const gr_float* FeS = spdens.contig1d_ptr(spid.dust(DustSpLUT::FeS_dust));
  const gr_float* Al2O3 = spdens.contig1d_ptr(spid.dust(DustSpLUT::Al2O3_dust));
  const gr_float* reforg = spdens.contig1d_ptr(spid.dust(DustSpLUT::ref_org_dust));
  const gr_float* volorg = spdens.contig1d_ptr(spid.dust(DustSpLUT::vol_org_dust));
  const gr_float* H2Oice = spdens.contig1d_ptr(spid.dust(DustSpLUT::H2O_ice_dust));

  // locals
  const double* const* kcol_buf = FullRxnRateBuf_kcol_bufs(&rxn_rate_buf);
  const double* const* kph_buf = FullRxnRateBuf_kph_bufs(&rxn_rate_buf);
  const double* const* grain_growth_rates =
    FullRxnRateBuf_grain_growth_bufs(&rxn_rate_buf);
  const double* h2dust = FullRxnRateBuf_h2dust(&rxn_rate_buf);
  // flag for if Gen Chiaki's dust model is enabled
  // -> to simplify bookkeeping in other parts of the code, we ignore whether
  //    the grain_growth flag was enabled, (the logic for setting the relevant
  //    reaction rates, explicitly sets the rates to 0 if disabled)
  const bool chiaki_model_dust_evolution = my_chemistry->dust_chemistry == 2;

  int i;
  double scoef, acoef;

  const int j = idx_range.j;
  const int k = idx_range.k;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI

        scoef  = kcol_buf[CollisionalRxnLUT::k2][i]*HII[i]*de[i];
        acoef  = kcol_buf[CollisionalRxnLUT::k1][i]*de[i]
               + kcol_buf[CollisionalRxnLUT::k57][i]*HI[i]
               + kcol_buf[CollisionalRxnLUT::k58][i]*HeI[i]/4.
               + kph_buf[PhotoRxnLUT::k24][i];
        out_spdens(spid.prim(PrimordialSpLUT::HI), i) =
            (scoef*dtit[i] + HI[i]) /
            (1. + acoef*dtit[i]);
        if (out_spdens(spid.prim(PrimordialSpLUT::HI), i) != out_spdens(spid.prim(PrimordialSpLUT::HI), i))  {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens(spid.prim(PrimordialSpLUT::HI), i),
                   HI[i],
                   HII[i],
                   de[i],
                   scoef,
                   acoef,
                   dtit [ i ]);
          }
          // ERROR_MESSAGE
        }

        // 2) HII
        scoef  = kcol_buf[CollisionalRxnLUT::k1][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*de[i]
               + kcol_buf[CollisionalRxnLUT::k57][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*out_spdens(spid.prim(PrimordialSpLUT::HI), i)
               + kcol_buf[CollisionalRxnLUT::k58][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*HeI[i]/4.
               + kph_buf[PhotoRxnLUT::k24][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i);
        acoef  = kcol_buf[CollisionalRxnLUT::k2][i]*de[i];
        out_spdens(spid.prim(PrimordialSpLUT::HII), i) =
            (scoef*dtit[i] + HII[i]) /
            (1. +acoef*dtit[i]);
        // 
        if (out_spdens(spid.prim(PrimordialSpLUT::HII), i) <= 0.)   //#####
        {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("negative HIIp! ::  %d %d %d %g %g %g %g %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens(spid.prim(PrimordialSpLUT::HII), i),
                   scoef,
                   dtit [ i ],
                   HII[i],
                   acoef,
                   kcol_buf[CollisionalRxnLUT::k2] [ i ],
                   de[i],
                   out_spdens(spid.prim(PrimordialSpLUT::HI), i),
                   kph_buf[PhotoRxnLUT::k24] [ i ]);
          }
        }

        // 3) Electron density

        scoef = 0.
                   + kcol_buf[CollisionalRxnLUT::k57][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*out_spdens(spid.prim(PrimordialSpLUT::HI), i)
                   + kcol_buf[CollisionalRxnLUT::k58][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*HeI[i]/4.
                   + kph_buf[PhotoRxnLUT::k24][i]*HI[i]
                   + kph_buf[PhotoRxnLUT::k25][i]*HeII[i]/4.
                   + kph_buf[PhotoRxnLUT::k26][i]*HeI[i]/4.;




        acoef = -(kcol_buf[CollisionalRxnLUT::k1][i]*HI[i]      - kcol_buf[CollisionalRxnLUT::k2][i]*HII[i]
                + kcol_buf[CollisionalRxnLUT::k3][i]*HeI[i]/4. -
             kcol_buf[CollisionalRxnLUT::k6][i]*HeIII[i]/4.
                + kcol_buf[CollisionalRxnLUT::k5][i]*HeII[i]/4. -
             kcol_buf[CollisionalRxnLUT::k4][i]*HeII[i]/4.);
        out_spdens(spid.prim(PrimordialSpLUT::e), i) =
            (scoef*dtit[i] + de[i]) /
            (1. + acoef*dtit[i]);

      }
    }

  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---

  for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE)  {

      // 4) HeI

      scoef  = kcol_buf[CollisionalRxnLUT::k4][i]*HeII[i]*de[i];
      acoef  = kcol_buf[CollisionalRxnLUT::k3][i]*de[i]
                   + kph_buf[PhotoRxnLUT::k26][i];

      if (my_chemistry->primordial_chemistry > 3)  {
        scoef = scoef +  4. * ( 0.
            + kcol_buf[CollisionalRxnLUT::k152][i] * HeHII[i] *    HI[i] /  5.
            + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII[i] *    de[i] /  5.
            );
        acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k148][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::k149][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::k150][i] *  H2II[i] /  2.;
      }
      out_spdens(spid.prim(PrimordialSpLUT::HeI), i) =
          ( scoef*dtit[i] + HeI[i] ) /
          ( 1. + acoef*dtit[i] );

      // 5) HeII

      scoef  = kcol_buf[CollisionalRxnLUT::k3][i]*out_spdens(spid.prim(PrimordialSpLUT::HeI), i)*de[i]
             + kcol_buf[CollisionalRxnLUT::k6][i]*HeIII[i]*de[i]
             + kph_buf[PhotoRxnLUT::k26][i]*out_spdens(spid.prim(PrimordialSpLUT::HeI), i);

      acoef  = kcol_buf[CollisionalRxnLUT::k4][i]*de[i] + kcol_buf[CollisionalRxnLUT::k5][i]*de[i]
             + kph_buf[PhotoRxnLUT::k25][i];
     
      if (my_chemistry->primordial_chemistry > 3)  {
        acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k151][i] *    HI[i];
      }
      out_spdens(spid.prim(PrimordialSpLUT::HeII), i) =
          ( scoef*dtit[i] + HeII[i] ) /
          ( 1. + acoef*dtit[i] );

      // 6) HeIII

      scoef   = kcol_buf[CollisionalRxnLUT::k5][i]*out_spdens(spid.prim(PrimordialSpLUT::HeII), i)*de[i]
              + kph_buf[PhotoRxnLUT::k25][i]*out_spdens(spid.prim(PrimordialSpLUT::HeII), i);
      acoef   = kcol_buf[CollisionalRxnLUT::k6][i]*de[i];
      out_spdens(spid.prim(PrimordialSpLUT::HeIII), i) =
          ( scoef*dtit[i] + HeIII[i] ) /
          ( 1. + acoef*dtit[i] );

    }
  }

  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms

    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {

        // 1) HI
        scoef  =      kcol_buf[CollisionalRxnLUT::k2][i] * HII[i] * de[i]
               + 2.*kcol_buf[CollisionalRxnLUT::k13][i]* HI[i]  * H2I[i]/2.
               +      kcol_buf[CollisionalRxnLUT::k11][i]* HII[i] * H2I[i]/2.
               + 2.*kcol_buf[CollisionalRxnLUT::k12][i]* de[i]  * H2I[i]/2.
               +      kcol_buf[CollisionalRxnLUT::k14][i]* HM[i]  * de[i]
               +      kcol_buf[CollisionalRxnLUT::k15][i]* HM[i]  * HI[i]
               + 2.*kcol_buf[CollisionalRxnLUT::k16][i]* HM[i]  * HII[i]
               + 2.*kcol_buf[CollisionalRxnLUT::k18][i]* H2II[i]* de[i]/2.
               +      kcol_buf[CollisionalRxnLUT::k19][i]* H2II[i]* HM[i]/2.
               + 2.*kph_buf[PhotoRxnLUT::k31][i]   * H2I[i]/2.;

        acoef  =      kcol_buf[CollisionalRxnLUT::k1][i] * de[i]
               +      kcol_buf[CollisionalRxnLUT::k7][i] * de[i]
               +      kcol_buf[CollisionalRxnLUT::k8][i] * HM[i]
               +      kcol_buf[CollisionalRxnLUT::k9][i] * HII[i]
               +      kcol_buf[CollisionalRxnLUT::k10][i]* H2II[i]/2.
               + 2.*kcol_buf[CollisionalRxnLUT::k22][i]* std::pow(HI[i],2)
               +      kcol_buf[CollisionalRxnLUT::k57][i]* HI[i]
               +      kcol_buf[CollisionalRxnLUT::k58][i]* HeI[i]/4.
               + kph_buf[PhotoRxnLUT::k24][i];

        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
            + kph_buf[PhotoRxnLUT::kdissHDI][i] * HDI[i]/3.0;
        }
        if ( (my_chemistry->metal_chemistry == 1)  &&
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
            + kph_buf[PhotoRxnLUT::kdissOH][i] * OHI[i] /17.0
            + kph_buf[PhotoRxnLUT::kdissH2O][i] * H2OI[i]/18.0;
        }

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i] != MASK_FALSE)  {
            acoef = acoef + 2. * h2dust[i] * rhoH[i];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf[CollisionalRxnLUT::k50][i] * HII[i] * DI[i]  / 2.
                + kcol_buf[CollisionalRxnLUT::k54][i] * H2I[i] * DI[i]  / 4.;
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k51][i] * DII[i] / 2.
                + kcol_buf[CollisionalRxnLUT::k55][i] * HDI[i] / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII[i] *    de[i] /  3.
              + kcol_buf[CollisionalRxnLUT::k134][i] *   HII[i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k135][i] *    HM[i] *    DI[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k150][i] *   HeI[i] *  H2II[i] /  8.
              + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII[i] *    de[i] /  5.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII[i] /  3.
              + kcol_buf[CollisionalRxnLUT::k130][i] *   DII[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k151][i] *  HeII[i] /  4.
              + kcol_buf[CollisionalRxnLUT::k152][i] * HeHII[i] /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI[i] *   H2I[i] / 24.
              + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI[i] *   H2I[i] / 32.
              + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII[i] *    OI[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I[i] *    CHI[i] / 26.
              + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I[i] *    OHI[i] / 34.
              + kcol_buf[CollisionalRxnLUT::kz26][i] *    OHI[i] *    COI[i] / 476.
              + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI[i] *    OHI[i] / 204.
              + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI[i] *    CHI[i] / 208.
              + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI[i] *    OHI[i] / 272.
              + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII[i] *    OHI[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII[i] *   H2OI[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII[i] *    O2I[i] / 32.
              + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII[i] *    OHI[i] / 204.
              + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII[i] *   H2I[i] / 32.
              + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII[i] *   H2I[i] / 34.
              + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII[i] *   H2I[i] / 36.
              + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII[i] *    de[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII[i] *    de[i] / 19.
              + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII[i] *    de[i] / 9.5
              + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI[i] *    OHI[i] / 476.
              + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI[i] *    OHI[i] / 748.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz15][i] *    CHI[i] / 13.
              + kcol_buf[CollisionalRxnLUT::kz16][i] *   CH2I[i] / 14.
              + kcol_buf[CollisionalRxnLUT::kz17][i] *    OHI[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz18][i] *   H2OI[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz19][i] *    O2I[i] / 32.
              + kcol_buf[CollisionalRxnLUT::kz27][i] *    CI[i] / 12.
              + kcol_buf[CollisionalRxnLUT::kz30][i] *    OI[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII[i] / 28.;
        }
        out_spdens(spid.prim(PrimordialSpLUT::HI), i) =
            ( scoef*dtit[i] + HI[i] ) /
            ( 1.f + acoef*dtit[i] );
        if (out_spdens(spid.prim(PrimordialSpLUT::HI), i) != out_spdens(spid.prim(PrimordialSpLUT::HI), i))  {
          OMP_PRAGMA_CRITICAL
          {
            std::printf("HUGE HIp! ::  %d %d %d %g %g %g %g %g\n",
                   i,
                   j,
                   k,
                   out_spdens(spid.prim(PrimordialSpLUT::HI), i),
                   HI[i],
                   HII[i],
                   de[i],
                   H2I[i]);
          }
        }

        // 2) HII

        scoef  =    kcol_buf[CollisionalRxnLUT::k1][i]  * HI[i] * de[i]
               +    kcol_buf[CollisionalRxnLUT::k10][i] * H2II[i]*HI[i]/2.
               +    kcol_buf[CollisionalRxnLUT::k57][i] * HI[i] * HI[i]
               +    kcol_buf[CollisionalRxnLUT::k58][i] * HI[i] * HeI[i]/4.
               + kph_buf[PhotoRxnLUT::k24][i]*HI[i];

        acoef  =    kcol_buf[CollisionalRxnLUT::k2][i]  * de[i]
               +    kcol_buf[CollisionalRxnLUT::k9][i]  * HI[i]
               +    kcol_buf[CollisionalRxnLUT::k11][i] * H2I[i]/2.
               +    kcol_buf[CollisionalRxnLUT::k16][i] * HM[i]
               +    kcol_buf[CollisionalRxnLUT::k17][i] * HM[i];
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf[CollisionalRxnLUT::k51][i] * HI[i] * DII[i] / 2.
                + kcol_buf[CollisionalRxnLUT::k52][i] * H2I[i] * DII[i] / 4.;
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k50][i] * DI[i] / 2.
                + kcol_buf[CollisionalRxnLUT::k53][i] * HDI[i] / 3.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII[i] *    HI[i] /  3.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k129][i] *    DI[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k134][i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k148][i] *   HeI[i] /  4.
              + kcol_buf[CollisionalRxnLUT::k149][i] *   HeI[i] /  4.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII[i] *    HI[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII[i] *    HI[i] / 28.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz22][i] *    OI[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz34][i] *    OHI[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz35][i] *   H2OI[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz36][i] *    O2I[i] / 32.;
        }
        out_spdens(spid.prim(PrimordialSpLUT::HII), i) =
            ( scoef*dtit[i] + HII[i] ) /
            ( 1. + acoef*dtit[i] );
        
        // 3) electrons:

        scoef =   kcol_buf[CollisionalRxnLUT::k8][i] * HM[i] * HI[i]
               +  kcol_buf[CollisionalRxnLUT::k15][i]* HM[i] * HI[i]
               +  kcol_buf[CollisionalRxnLUT::k17][i]* HM[i] * HII[i]
               +  kcol_buf[CollisionalRxnLUT::k57][i]* HI[i] * HI[i]
               +  kcol_buf[CollisionalRxnLUT::k58][i]* HI[i] * HeI[i]/4.
        // 
               + kph_buf[PhotoRxnLUT::k24][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)
               + kph_buf[PhotoRxnLUT::k25][i]*out_spdens(spid.prim(PrimordialSpLUT::HeII), i)/4.
               + kph_buf[PhotoRxnLUT::k26][i]*out_spdens(spid.prim(PrimordialSpLUT::HeI), i)/4.;

        if ( (my_chemistry->metal_chemistry == 1)  &&
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef
            + kph_buf[PhotoRxnLUT::kphCI][i] * CI[i]/12.0
            + kph_buf[PhotoRxnLUT::kphOI][i] * OI[i]/16.0;
        }

        acoef = - (kcol_buf[CollisionalRxnLUT::k1][i] *HI[i]    - kcol_buf[CollisionalRxnLUT::k2][i]*HII[i]
                +  kcol_buf[CollisionalRxnLUT::k3][i] *HeI[i]/4. -
             kcol_buf[CollisionalRxnLUT::k6][i]*HeIII[i]/4.
                +  kcol_buf[CollisionalRxnLUT::k5][i] *HeII[i]/4. -
             kcol_buf[CollisionalRxnLUT::k4][i]*HeII[i]/4.
                +  kcol_buf[CollisionalRxnLUT::k14][i]*HM[i]
                -  kcol_buf[CollisionalRxnLUT::k7][i] *HI[i]
                -  kcol_buf[CollisionalRxnLUT::k18][i]*H2II[i]/2.);
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef
                + kcol_buf[CollisionalRxnLUT::k56][i] * DI[i] * HM[i] / 2.;
          acoef = acoef
                - kcol_buf[CollisionalRxnLUT::k1] [i] * DI[i] / 2.
                + kcol_buf[CollisionalRxnLUT::k2] [i] * DII[i] / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM[i] *    HI[i] /  2.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII[i] /  3.
              + kcol_buf[CollisionalRxnLUT::k132][i] *    DI[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k153][i] * HeHII[i] /  5.;
        }

        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz44][i] *   CII[i] / 12.
              + kcol_buf[CollisionalRxnLUT::kz45][i] *   OII[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII[i] / 19.
              + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII[i] / 19.
              + kcol_buf[CollisionalRxnLUT::kz50][i] *  O2II[i] / 32.;
        }
        out_spdens(spid.prim(PrimordialSpLUT::e), i) =
            ( scoef*dtit[i] + de[i] ) /
            ( 1. + acoef*dtit[i] );

        // 7) H2

        scoef = 2.*(kcol_buf[CollisionalRxnLUT::k8][i]  * HM[i]   * HI[i]
              +       kcol_buf[CollisionalRxnLUT::k10][i] * H2II[i] * HI[i]/2.
              +       kcol_buf[CollisionalRxnLUT::k19][i] * H2II[i] * HM[i]/2.
              +       kcol_buf[CollisionalRxnLUT::k22][i] * HI[i] * std::pow((HI[i]),2.));
        acoef = ( kcol_buf[CollisionalRxnLUT::k13][i]*HI[i] + kcol_buf[CollisionalRxnLUT::k11][i]*HII[i]
                + kcol_buf[CollisionalRxnLUT::k12][i]*de[i] )
                + kph_buf[PhotoRxnLUT::k29][i] + kph_buf[PhotoRxnLUT::k31][i];

        if (anydust != MASK_FALSE)  {
          if(itmask_metal[i] != MASK_FALSE)  {
            scoef = scoef + 2. * h2dust[i] *
                 HI[i] * rhoH[i];
          }
        }
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          scoef = scoef + 2. * (
                  kcol_buf[CollisionalRxnLUT::k53][i] * HDI[i] * HII[i] / 3.
                + kcol_buf[CollisionalRxnLUT::k55][i] * HDI[i] * HI[i] / 3.
                   );
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k52][i] * DII[i] / 2.
                + kcol_buf[CollisionalRxnLUT::k54][i] * DI[i] / 2.;
        }
#endif
        if ( (my_chemistry->metal_chemistry == 1)  && 
             (itmask_metal[i] != MASK_FALSE) )  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI[i] *    CHI[i] / 13.
              + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI[i] *   CH2I[i] / 14.
              + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI[i] *    OHI[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI[i] *   H2OI[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII[i] *    de[i] / 18.
             );
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI[i] / 12.
              + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz23][i] *    CHI[i] / 13.
              + kcol_buf[CollisionalRxnLUT::kz24][i] *    OHI[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII[i] / 16.
              + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII[i] / 17.
              + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII[i] / 18.
              + kcol_buf[CollisionalRxnLUT::kz51][i] *    CI[i] / 12.;
          if (chiaki_model_dust_evolution)  {
            if (my_chemistry->dust_species > 0)  {
              scoef = scoef + 2. *
                    grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 2.;

            }
            if (my_chemistry->dust_species > 1)  {
              scoef = scoef + 2. * (
                    grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 3.
                  + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 4.
                  + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i]
                  + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] * 3.
                );
            }
            if (my_chemistry->dust_species > 2)  {
              acoef = acoef
              + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] / H2I[i] * 2. * 2.;
            }
          }
        }
        out_spdens(spid.prim(PrimordialSpLUT::H2I), i) =
            ( scoef*dtit[i] + H2I[i] ) /
            ( 1. + acoef*dtit[i] );

        // 8) H-

        scoef = kcol_buf[CollisionalRxnLUT::k7][i] * HI[i] * de[i];
        acoef = (kcol_buf[CollisionalRxnLUT::k8][i]  + kcol_buf[CollisionalRxnLUT::k15][i])  * HI[i] +
                (kcol_buf[CollisionalRxnLUT::k16][i] + kcol_buf[CollisionalRxnLUT::k17][i])  * HII[i] +
               kcol_buf[CollisionalRxnLUT::k14][i] * de[i] + kcol_buf[CollisionalRxnLUT::k19][i] * H2II[i]/2.0f +
               kph_buf[PhotoRxnLUT::k27][i];
#ifdef CONTRIBUTION_OF_MINOR_SPECIES
        if (my_chemistry->primordial_chemistry > 2)  {
          acoef = acoef
                + kcol_buf[CollisionalRxnLUT::k56][i] * DI[i] / 2.;
        }
#endif
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM[i] *    HI[i] /  2.;
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k135][i] *    DI[i] /  2.;
        }
        out_spdens(spid.prim(PrimordialSpLUT::HM), i) =
            (scoef*dtit[i] + HM[i]) /
            (1.0f + acoef*dtit[i]);


        // 9) H2+

        out_spdens(spid.prim(PrimordialSpLUT::H2II), i) =
            2. * (
                kcol_buf[CollisionalRxnLUT::k9][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
              + kcol_buf[CollisionalRxnLUT::k11][i]*out_spdens(spid.prim(PrimordialSpLUT::H2I), i)/2.*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
              + kcol_buf[CollisionalRxnLUT::k17][i]*out_spdens(spid.prim(PrimordialSpLUT::HM), i)*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
              + kph_buf[PhotoRxnLUT::k29][i]*out_spdens(spid.prim(PrimordialSpLUT::H2I), i)
                 )
            / (
                kcol_buf[CollisionalRxnLUT::k10][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)
              + kcol_buf[CollisionalRxnLUT::k18][i]*out_spdens(spid.prim(PrimordialSpLUT::e), i)
              + kcol_buf[CollisionalRxnLUT::k19][i]*out_spdens(spid.prim(PrimordialSpLUT::HM), i)
              + (kph_buf[PhotoRxnLUT::k28][i]+kph_buf[PhotoRxnLUT::k30][i])
              );
        if (my_chemistry->primordial_chemistry > 3)  {
          out_spdens(spid.prim(PrimordialSpLUT::H2II), i) =
              2. * (
                  kcol_buf[CollisionalRxnLUT::k9] [i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
                + kcol_buf[CollisionalRxnLUT::k11][i]*out_spdens(spid.prim(PrimordialSpLUT::H2I), i)/2.*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
                + kcol_buf[CollisionalRxnLUT::k17][i]*out_spdens(spid.prim(PrimordialSpLUT::HM), i)*out_spdens(spid.prim(PrimordialSpLUT::HII), i)
                + kph_buf[PhotoRxnLUT::k29][i]*out_spdens(spid.prim(PrimordialSpLUT::H2I), i)
                + kcol_buf[CollisionalRxnLUT::k152][i]*HeHII[i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)/5.
                   )
              / (
                  kcol_buf[CollisionalRxnLUT::k10][i]*out_spdens(spid.prim(PrimordialSpLUT::HI), i)
                + kcol_buf[CollisionalRxnLUT::k18][i]*out_spdens(spid.prim(PrimordialSpLUT::e), i)
                + kcol_buf[CollisionalRxnLUT::k19][i]*out_spdens(spid.prim(PrimordialSpLUT::HM), i)
                + (kph_buf[PhotoRxnLUT::k28][i]+kph_buf[PhotoRxnLUT::k30][i])
                + kcol_buf[CollisionalRxnLUT::k150][i]*out_spdens(spid.prim(PrimordialSpLUT::HeI), i)/4.
                );
        }
      }
    }
    // 
  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {
        
        // 1) DI
        scoef =   (       kcol_buf[CollisionalRxnLUT::k2][i] * DII[i] * de[i]
                   +      kcol_buf[CollisionalRxnLUT::k51][i]* DII[i] * HI[i]
                   + 2.*kcol_buf[CollisionalRxnLUT::k55][i]* HDI[i] * HI[i]/3.
                   )
                   + 2. * kph_buf[PhotoRxnLUT::kdissHDI][i] * HDI[i]/3.0;
        acoef  =    kcol_buf[CollisionalRxnLUT::k1][i] * de[i]
               +    kcol_buf[CollisionalRxnLUT::k50][i] * HII[i]
               +    kcol_buf[CollisionalRxnLUT::k54][i] * H2I[i]/2.
               +    kcol_buf[CollisionalRxnLUT::k56][i] * HM[i]
               + kph_buf[PhotoRxnLUT::k24][i];
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  2. * ( 0.
              + kcol_buf[CollisionalRxnLUT::k131][i] *  HDII[i] *    de[i] /  3.
              + kcol_buf[CollisionalRxnLUT::k133][i] *   DII[i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k134][i] *   HII[i] *    DM[i] /  2.
              + kcol_buf[CollisionalRxnLUT::k136][i] *    DM[i] *    HI[i] /  2.
              );
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k129][i] *   HII[i]
              + kcol_buf[CollisionalRxnLUT::k132][i] *    de[i]
              + kcol_buf[CollisionalRxnLUT::k135][i] *    HM[i];
        }
        out_spdens(spid.prim(PrimordialSpLUT::DI), i) =
            ( scoef*dtit[i] + DI[i] ) /
            ( 1. + acoef*dtit[i] );

        // 2) DII
        scoef =   (   kcol_buf[CollisionalRxnLUT::k1][i]  * DI[i] * de[i]
              +       kcol_buf[CollisionalRxnLUT::k50][i] * HII[i]* DI[i]
              +  2.*kcol_buf[CollisionalRxnLUT::k53][i] * HII[i]* HDI[i]/3.
              )
              + kph_buf[PhotoRxnLUT::k24][i]*DI[i];
        acoef = 0.;
        acoef =    kcol_buf[CollisionalRxnLUT::k2][i]  * de[i]
              +    kcol_buf[CollisionalRxnLUT::k51][i] * HI[i]
              +    kcol_buf[CollisionalRxnLUT::k52][i] * H2I[i]/2.;
        if (my_chemistry->primordial_chemistry > 3)  {
          acoef = acoef
              + kcol_buf[CollisionalRxnLUT::k130][i] *    HI[i]
              + kcol_buf[CollisionalRxnLUT::k133][i] *    DM[i] /  2.;
        }
        out_spdens(spid.prim(PrimordialSpLUT::DII), i) =
            ( scoef*dtit[i] + DII[i] ) /
            ( 1. + acoef*dtit[i] );

        // 3) HDI
        scoef = 3.*(kcol_buf[CollisionalRxnLUT::k52][i] * DII[i]*
             H2I[i]/2./2.
             + kcol_buf[CollisionalRxnLUT::k54][i] * DI[i] * H2I[i]/2./2.
        // !   &           + 2._DKIND*k56(i) * DI[i] * HM[i]/2._DKIND
        //- ! corrected by GC202005
             +          kcol_buf[CollisionalRxnLUT::k56][i] * DI[i] * HM[i]/2.
                   );
        acoef  =    kcol_buf[CollisionalRxnLUT::k53][i] * HII[i]
               +    kcol_buf[CollisionalRxnLUT::k55][i] * HI[i]
               + kph_buf[PhotoRxnLUT::kdissHDI][i];
        if (my_chemistry->primordial_chemistry > 3)  {
          scoef = scoef +  3. * ( 0.
              + kcol_buf[CollisionalRxnLUT::k125][i] *  HDII[i] *    HI[i] /  3.
              + kcol_buf[CollisionalRxnLUT::k137][i] *    DM[i] *    HI[i] /  2.
              );
        }
        out_spdens(spid.prim(PrimordialSpLUT::HDI), i) =
            ( scoef*dtit[i] + HDI[i] ) /
            ( 1. + acoef*dtit[i] );

      }
    }
  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask[i] != MASK_FALSE)  {
        
        // 1) DM
        scoef =
              kcol_buf[CollisionalRxnLUT::k132][i] *    DI[i] *    de[i]
            + kcol_buf[CollisionalRxnLUT::k135][i] *    HM[i] *    DI[i];
        acoef =
              kcol_buf[CollisionalRxnLUT::k133][i] *   DII[i] /  2.
            + kcol_buf[CollisionalRxnLUT::k134][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::k136][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::k137][i] *    HI[i];

        out_spdens(spid.prim(PrimordialSpLUT::DM), i) =
            ( scoef*dtit[i] + DM[i] ) /
            ( 1. + acoef*dtit[i] );

        // 2) HDII
        scoef = 3. * (
              kcol_buf[CollisionalRxnLUT::k129][i] *    DI[i] *   HII[i] /  2.
            + kcol_buf[CollisionalRxnLUT::k130][i] *   DII[i] *    HI[i] /  2.
            );
        acoef =
              kcol_buf[CollisionalRxnLUT::k125][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::k131][i] *    de[i];

        out_spdens(spid.prim(PrimordialSpLUT::HDII), i) =
            ( scoef*dtit[i] + HDII[i] ) /
            ( 1. + acoef*dtit[i] );

        // 3) HeHII
        scoef = 5. * (
              kcol_buf[CollisionalRxnLUT::k148][i] *   HeI[i] *   HII[i] /  4.
            + kcol_buf[CollisionalRxnLUT::k149][i] *   HeI[i] *   HII[i] /  4.
            + kcol_buf[CollisionalRxnLUT::k150][i] *   HeI[i] *  H2II[i] /  8.
            + kcol_buf[CollisionalRxnLUT::k151][i] *  HeII[i] *    HI[i] /  4.
            );
        acoef =
              kcol_buf[CollisionalRxnLUT::k152][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::k153][i] *    de[i];

        out_spdens(spid.prim(PrimordialSpLUT::HeHII), i) =
            ( scoef*dtit[i] + HeHII[i] ) /
            ( 1. + acoef*dtit[i] );

      }
    }
  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE)  {

        // ***** CI **********
        scoef = 0. + 12. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI[i] *    CHI[i] / 13.
            + kcol_buf[CollisionalRxnLUT::kz44][i] *   CII[i] *    de[i] / 12.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz20][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz27][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    OHI[i] / 17.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    O2I[i] / 32.
            + kcol_buf[CollisionalRxnLUT::kz51][i] *   H2I[i] /  2.
            + kph_buf[PhotoRxnLUT::kphCI][i];
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::AC_dust]      [i] / CI[i] * 12.;
          }
        }
        scoef = scoef + 12. *
            kph_buf[PhotoRxnLUT::kdissCO][i] * COI[i] /28.0;

        out_spdens(spid.metal(MetalSpLUT::CI), i) =
            ( scoef*dtit[i] + CI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** CII **********
        scoef = 0. + 12. * ( 0.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz37][i] *    OHI[i] / 17.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *    O2I[i] / 32.
            + kcol_buf[CollisionalRxnLUT::kz44][i] *    de[i];
        scoef = scoef
            + kph_buf[PhotoRxnLUT::kphCI][i] * CI[i];

        out_spdens(spid.metal(MetalSpLUT::CII), i) =
            ( scoef*dtit[i] + CII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** COI *********
        scoef = 0. + 28. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI[i] *    OHI[i] / 204.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI[i] *    O2I[i] / 384.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI[i] *    CHI[i] / 208.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII[i] *    O2I[i] / 384.
            + kcol_buf[CollisionalRxnLUT::kz43][i] *  COII[i] *    HI[i] / 28.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    OHI[i] / 17.;
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] / COI[i] * 17. * 0.5
            + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] / COI[i] * 17.;
          }
        }
        acoef = acoef
            + kph_buf[PhotoRxnLUT::kdissCO][i];

        out_spdens(spid.metal(MetalSpLUT::COI), i) =
            ( scoef*dtit[i] + COI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** CO2I *********
        scoef = 0. + 44. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    OHI[i] *    COI[i] / 476.
           );
        acoef = 0.;

        out_spdens(spid.metal(MetalSpLUT::CO2I), i) =
            ( scoef*dtit[i] + CO2I[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** OI **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI[i] *    OHI[i] / 17.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI[i] *    O2I[i] / 32.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OHI[i] *    OHI[i] / 289.
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI[i] *    O2I[i] / 384.
            + kcol_buf[CollisionalRxnLUT::kz39][i] *   OII[i] *    HI[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz45][i] *   OII[i] *    de[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz47][i] * H2OII[i] *    de[i] / 18.
            + kcol_buf[CollisionalRxnLUT::kz50][i] *  O2II[i] *    de[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI[i] *    O2I[i] / 896.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz21][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::kz30][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz31][i] *    OI[i] / 8.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    CHI[i] / 13.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OHI[i] / 17.
            + kph_buf[PhotoRxnLUT::kphOI][i];
        scoef = scoef + 16. *
            ( kph_buf[PhotoRxnLUT::kdissOH][i] * OHI[i] /17.0
            + kph_buf[PhotoRxnLUT::kdissCO][i] * COI[i] /28.0);

        out_spdens(spid.metal(MetalSpLUT::OI), i) =
            ( scoef*dtit[i] + OI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** OH *********
        scoef = 0. + 17. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI[i] *   H2OI[i] / 18.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI[i] *    O2I[i] / 32.
            + kcol_buf[CollisionalRxnLUT::kz21][i] *    OI[i] *   H2I[i] / 32.
            + kcol_buf[CollisionalRxnLUT::kz30][i] *    OI[i] *    HI[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz46][i] * H2OII[i] *    de[i] / 18.
            + kcol_buf[CollisionalRxnLUT::kz49][i] * H3OII[i] *    de[i] / 19.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz17][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OHI[i] / 8.5
            + kcol_buf[CollisionalRxnLUT::kz26][i] *    COI[i] / 28.
            + kcol_buf[CollisionalRxnLUT::kz28][i] *    CI[i] / 12.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII[i] / 12.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI[i] / 28.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI[i] / 44.;
        acoef = acoef
            + kph_buf[PhotoRxnLUT::kdissOH][i];
        scoef = scoef + 17. *
            kph_buf[PhotoRxnLUT::kdissH2O][i] * H2OI[i]/18.0;

        out_spdens(spid.metal(MetalSpLUT::OHI), i) =
            ( scoef*dtit[i] + OHI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** H2O *********
        scoef = 0. + 18. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz24][i] *   H2I[i] *    OHI[i] / 34.
            + kcol_buf[CollisionalRxnLUT::kz25][i] *    OHI[i] *    OHI[i] / 289.
            + kcol_buf[CollisionalRxnLUT::kz48][i] * H3OII[i] *    de[i] / 19.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz18][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII[i];
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / H2OI[i] * 18. * 2.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / H2OI[i] * 18. * 3.
            + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] / H2OI[i] * 18. * 4.
            + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] / H2OI[i] * 18.
            + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] / H2OI[i] * 18. * 3.;
          }
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust]  [i] / H2OI[i] * 18.;
          }
        }
        acoef = acoef
            + kph_buf[PhotoRxnLUT::kdissH2O][i];

        out_spdens(spid.metal(MetalSpLUT::H2OI), i) =
            ( scoef*dtit[i] + H2OI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** O2I *********
        scoef = 0. + 32. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz31][i] *    OI[i] *    OI[i] / 256.
            + kcol_buf[CollisionalRxnLUT::kz33][i] *    OI[i] *    OHI[i] / 272.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz19][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz29][i] *    CI[i] / 12.
            + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII[i]
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII[i] / 12.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI[i] / 28.;

        out_spdens(spid.metal(MetalSpLUT::O2I), i) =
            ( scoef*dtit[i] + O2I[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** SiI **********
        scoef = 0. + 28. * ( 0.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *    OHI[i] / 17.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *    O2I[i] / 32.;
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::SiM_dust]     [i] / SiI[i] * 28.;
          }
        }

        out_spdens(spid.metal(MetalSpLUT::SiI), i) =
            ( scoef*dtit[i] + SiI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** SiOI **********
        scoef = 0. + 44. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz52][i] *   SiI[i] *    OHI[i] / 476.
            + kcol_buf[CollisionalRxnLUT::kz53][i] *   SiI[i] *    O2I[i] / 896.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *    OHI[i] / 17.;
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 0)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / SiOI[i] * 44.;
          }
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / SiOI[i] * 44.;
          }
        }

        out_spdens(spid.metal(MetalSpLUT::SiOI), i) =
            ( scoef*dtit[i] + SiOI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** SiO2I **********
        scoef = 0. + 60. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz54][i] *  SiOI[i] *    OHI[i] / 748.
           );
        acoef = 0.;
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust]   [i] / SiO2I[i] * 60.;
          }
        }

        out_spdens(spid.metal(MetalSpLUT::SiO2I), i) =
            ( scoef*dtit[i] + SiO2I[i] ) /
            ( 1. + acoef*dtit[i] );

        // MINOR BUT IMPORTANT SPECIES FOR MOLECULAR FORMATION
        //- ***** CH **********
        scoef = 0. + 13. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI[i] *   CH2I[i] / 14.
            + kcol_buf[CollisionalRxnLUT::kz20][i] *    CI[i] *   H2I[i] / 24.
            + kcol_buf[CollisionalRxnLUT::kz27][i] *    CI[i] *    HI[i] / 12.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz15][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz32][i] *    OI[i] / 16.;

        out_spdens(spid.metal(MetalSpLUT::CHI), i) =
            ( scoef*dtit[i] + CHI[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** CH2I *********
        scoef = 0. + 14. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz23][i] *   H2I[i] *    CHI[i] / 26.
            + kcol_buf[CollisionalRxnLUT::kz51][i] *   H2I[i] *    CI[i] / 24.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz16][i] *    HI[i];
        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 2)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] / CH2I[i] * 14. * 0.5;
          }
        }

        out_spdens(spid.metal(MetalSpLUT::CH2I), i) =
            ( scoef*dtit[i] + CH2I[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** COII **********
        scoef = 0. + 28. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz37][i] *   CII[i] *    OHI[i] / 204.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz43][i] *    HI[i];

        out_spdens(spid.metal(MetalSpLUT::COII), i) =
            ( scoef*dtit[i] + COII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** OII **********
        scoef = 0. + 16. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz22][i] *   HII[i] *    OI[i] / 16.
            + kcol_buf[CollisionalRxnLUT::kz38][i] *   CII[i] *    O2I[i] / 384.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz39][i] *    HI[i]
            + kcol_buf[CollisionalRxnLUT::kz40][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz45][i] *    de[i];
        scoef = scoef
            + kph_buf[PhotoRxnLUT::kphOI][i] * OI[i];

        out_spdens(spid.metal(MetalSpLUT::OII), i) =
            ( scoef*dtit[i] + OII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** OHII **********
        scoef = 0. + 17. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz34][i] *   HII[i] *    OHI[i] / 17.
            + kcol_buf[CollisionalRxnLUT::kz40][i] *   OII[i] *   H2I[i] / 32.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz41][i] *   H2I[i] /  2.;

        out_spdens(spid.metal(MetalSpLUT::OHII), i) =
            ( scoef*dtit[i] + OHII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** H2OII **********
        scoef = 0. + 18. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz35][i] *   HII[i] *   H2OI[i] / 18.
            + kcol_buf[CollisionalRxnLUT::kz41][i] *  OHII[i] *   H2I[i] / 34.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz42][i] *   H2I[i] /  2.
            + kcol_buf[CollisionalRxnLUT::kz46][i] *    de[i]
            + kcol_buf[CollisionalRxnLUT::kz47][i] *    de[i];

        out_spdens(spid.metal(MetalSpLUT::H2OII), i) =
            ( scoef*dtit[i] + H2OII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** H3OII **********
        scoef = 0. + 19. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz42][i] * H2OII[i] *   H2I[i] / 36.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz48][i] *    de[i]
            + kcol_buf[CollisionalRxnLUT::kz49][i] *    de[i];

        out_spdens(spid.metal(MetalSpLUT::H3OII), i) =
            ( scoef*dtit[i] + H3OII[i] ) /
            ( 1. + acoef*dtit[i] );


        // ***** O2II **********
        scoef = 0. + 32. * ( 0.
            + kcol_buf[CollisionalRxnLUT::kz36][i] *   HII[i] *    O2I[i] / 32.
           );
        acoef = 0.
            + kcol_buf[CollisionalRxnLUT::kz50][i] *    de[i];

        out_spdens(spid.metal(MetalSpLUT::O2II), i) =
            ( scoef*dtit[i] + O2II[i] ) /
            ( 1. + acoef*dtit[i] );


        if (chiaki_model_dust_evolution)  {
          if (my_chemistry->dust_species > 0)  {
            // ***** Mg **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] / Mg[i] * 24.;
            if (my_chemistry->dust_species > 1)  {
              acoef = acoef
              + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] / Mg[i] * 24. * 2.
              + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] / Mg[i] * 24.;
            }

            out_spdens(spid.metal(MetalSpLUT::Mg), i) =
                ( scoef*dtit[i] + Mg[i] ) /
                ( 1. + acoef*dtit[i] );

          }

          if (my_chemistry->dust_species > 1)  {
            // ***** Al **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] / Al[i] * 27. * 2.;

            out_spdens(spid.metal(MetalSpLUT::Al), i) =
                ( scoef*dtit[i] + Al[i] ) /
                ( 1. + acoef*dtit[i] );


            // ***** S  **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] / S[i] * 32.;

            out_spdens(spid.metal(MetalSpLUT::S), i) =
                ( scoef*dtit[i] + S[i] ) /
                ( 1. + acoef*dtit[i] );


            // ***** Fe **********
            scoef = 0.;
            acoef = 0.;
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::FeM_dust]     [i] / Fe[i] * 56.
            + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] / Fe[i] * 56. * 3.
            + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] / Fe[i] * 56.;

            out_spdens(spid.metal(MetalSpLUT::Fe), i) =
                ( scoef*dtit[i] + Fe[i] ) /
                ( 1. + acoef*dtit[i] );

          }
        }

      }
    }
  }

  // --- (D4) Now do dust species ---
  if (chiaki_model_dust_evolution)  {
    for (i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE)  {

        if (my_chemistry->dust_species > 0)  {
          // ***** MgSiO3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust]  [i] * 100.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::MgSiO3_dust), i) =
              ( scoef*dtit[i] + MgSiO3[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** AC **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::AC_dust]      [i] * 12.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::AC_dust), i) =
              ( scoef*dtit[i] + AC[i] ) /
              ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 1)  {
          // ***** SiM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::SiM_dust]     [i] * 28.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::SiM_dust), i) =
              ( scoef*dtit[i] + SiM[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** FeM **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::FeM_dust]     [i] * 56.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::FeM_dust), i) =
              ( scoef*dtit[i] + FeM[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** Mg2SiO4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust] [i] * 140.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::Mg2SiO4_dust), i) =
              ( scoef*dtit[i] + Mg2SiO4[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** Fe3O4 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust]   [i] * 232.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::Fe3O4_dust), i) =
              ( scoef*dtit[i] + Fe3O4[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** SiO2D **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust]   [i] * 60.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::SiO2_dust), i) =
              ( scoef*dtit[i] + SiO2D[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** MgO **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::MgO_dust]     [i] * 40.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::MgO_dust), i) =
              ( scoef*dtit[i] + MgO[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** FeS **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust]     [i] * 88.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::FeS_dust), i) =
              ( scoef*dtit[i] + FeS[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** Al2O3 **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust]   [i] * 102.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::Al2O3_dust), i) =
              ( scoef*dtit[i] + Al2O3[i] ) /
              ( 1. + acoef*dtit[i] );

        }

        if (my_chemistry->dust_species > 2)  {
          // ***** reforg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust]  [i] * 22.68;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::ref_org_dust), i) =
              ( scoef*dtit[i] + reforg[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** volorg **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust]  [i] * 32.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::vol_org_dust), i) =
              ( scoef*dtit[i] + volorg[i] ) /
              ( 1. + acoef*dtit[i] );


          // ***** H2Oice **********
          scoef = 0.;
          scoef = scoef
          + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust]  [i] * 18.;
          acoef = 0.;

          out_spdens(spid.dust(DustSpLUT::H2O_ice_dust), i) =
              ( scoef*dtit[i] + H2Oice[i] ) /
              ( 1. + acoef*dtit[i] );

        }

      }
    }
  }
}


/// Computes the time derivatives of the mass densities of various species
///
/// This function is specialized to operate on a single cell at a time. This is
/// mostly for historical reasons. In reality, it would be faster to operate on
/// many species at once.
///
/// @param[out] deriv This is where we store the time derivative of the
///     species mass densities. The caller must make sure this is filled with
///     zeros before calling this function
/// @param[in]  anydust
/// @param[in]  rhoH Indicates the mass density of all Hydrogen
/// @param[in]  itmask_metal Indicates where we should account for metal
///     chemistry
/// @param[in]  my_chemistry Provides various runtime parameters (we probably
///     don't need to pass the whole thing)
/// @param[in]  spdens Specifies the densities of the various species
///     that Grackle evolves in a format that allows the values to be accessed
///     with the index lookup table.
/// @param[in]  rxn_rate_buf specifies the precomputed rxn rates (depends on
///     local physical conditions)
/// @param[in]  sp_kind_map associates intervals of species ids with the kind
///     with the appropriate @ref SpKind
///
/// @note
/// Some of occurences of the const keyword is somewhat symbolic (to indicate
/// that we won't modify an arg or the contents of a struct). In reality,
/// C/C++ allows you to modify pointers in a `const struct`
///
/// @todo
/// We should consider how to deduplicate this logic and most of the logic in
/// `step_rate_g`
/// -> I am pretty sure Gen performed a copy-and-paste when he wrote this
/// -> it is totally unsustainable for us to maintain both implementations
/// -> A common obstacle for unifying this logic is:
///    - step_rate_g uses these calculations with a backward difference
///      formula (see 3.2 of the grackle paper)
///    - Here, we're computing the net time derivative.
///    We can DEFINITELY overcome this obstacle!
/// -> there's another obstacle that applies for unifying another part of
///    the logic:
///    - Parts (A), (B), and (C) in the original `step_rate_g` all apply
///      partial updates (again, see 3.2 of the grackle paper)
///    - notably, part (D) does not do this (but I'm not sure that there
///      is a reasonable motivation for that choice)
inline void species_density_derivatives_0d(
  View<gr_float**> deriv, gr_mask_type anydust,
  const double* rhoH, const gr_mask_type* itmask_metal,
  const chemistry_data* my_chemistry,
  SpeciesMultiView<const gr_float> spdens,
  const FullRxnRateBuf rxn_rate_buf, const PartMap<SpKind>& sp_kind_map
) {
  const SpIdAdaptor spid(sp_kind_map);

  // define some local variables carried over from the fortran version:
  // - the goal is to eventually remove all of these
  // - originally we only conditionally defined some of these variables, but
  //   there honestly isn't any benefit to doing that (the memory is allocated
  //   already)
  // - originally, each of these variables was a stack allocated variable that
  //   held a copy of the corresponding entry in dsp. Now, these variables are
  //   references to casted copies taken from dsp
  // - we aren't being that consistent with the historic implementation when
  //   gr_float isn't the same as double (but, the historical implementation
  //   definitely didn't handle this case properly)

  const gr_float& de      = spdens(spid.prim(PrimordialSpLUT::e), 0);
  const gr_float& HI      = spdens(spid.prim(PrimordialSpLUT::HI), 0);
  const gr_float& HII     = spdens(spid.prim(PrimordialSpLUT::HII), 0);
  const gr_float& HeI     = spdens(spid.prim(PrimordialSpLUT::HeI), 0);
  const gr_float& HeII    = spdens(spid.prim(PrimordialSpLUT::HeII), 0);
  const gr_float& HeIII   = spdens(spid.prim(PrimordialSpLUT::HeIII), 0);
  const gr_float& HM      = spdens(spid.prim(PrimordialSpLUT::HM), 0);
  const gr_float& H2I     = spdens(spid.prim(PrimordialSpLUT::H2I), 0);
  const gr_float& H2II    = spdens(spid.prim(PrimordialSpLUT::H2II), 0);
  const gr_float& DI      = spdens(spid.prim(PrimordialSpLUT::DI), 0);
  const gr_float& DII     = spdens(spid.prim(PrimordialSpLUT::DII), 0);
  const gr_float& HDI     = spdens(spid.prim(PrimordialSpLUT::HDI), 0);
  const gr_float& DM      = spdens(spid.prim(PrimordialSpLUT::DM), 0);
  const gr_float& HDII    = spdens(spid.prim(PrimordialSpLUT::HDII), 0);
  const gr_float& HeHII   = spdens(spid.prim(PrimordialSpLUT::HeHII), 0);
  const gr_float& CI      = spdens(spid.metal(MetalSpLUT::CI), 0);
  const gr_float& CII     = spdens(spid.metal(MetalSpLUT::CII), 0);
  const gr_float& COI     = spdens(spid.metal(MetalSpLUT::COI), 0);
  const gr_float& CO2I    = spdens(spid.metal(MetalSpLUT::CO2I), 0);
  const gr_float& OI      = spdens(spid.metal(MetalSpLUT::OI), 0);
  const gr_float& OHI     = spdens(spid.metal(MetalSpLUT::OHI), 0);
  const gr_float& H2OI    = spdens(spid.metal(MetalSpLUT::H2OI), 0);
  const gr_float& O2I     = spdens(spid.metal(MetalSpLUT::O2I), 0);
  const gr_float& SiI     = spdens(spid.metal(MetalSpLUT::SiI), 0);
  const gr_float& SiOI    = spdens(spid.metal(MetalSpLUT::SiOI), 0);
  const gr_float& SiO2I   = spdens(spid.metal(MetalSpLUT::SiO2I), 0);
  const gr_float& CHI     = spdens(spid.metal(MetalSpLUT::CHI), 0);
  const gr_float& CH2I    = spdens(spid.metal(MetalSpLUT::CH2I), 0);
  const gr_float& COII    = spdens(spid.metal(MetalSpLUT::COII), 0);
  const gr_float& OII     = spdens(spid.metal(MetalSpLUT::OII), 0);
  const gr_float& OHII    = spdens(spid.metal(MetalSpLUT::OHII), 0);
  const gr_float& H2OII   = spdens(spid.metal(MetalSpLUT::H2OII), 0);
  const gr_float& H3OII   = spdens(spid.metal(MetalSpLUT::H3OII), 0);
  const gr_float& O2II    = spdens(spid.metal(MetalSpLUT::O2II), 0);
  const gr_float& Mg      = spdens(spid.metal(MetalSpLUT::Mg), 0);
  const gr_float& Al      = spdens(spid.metal(MetalSpLUT::Al), 0);
  const gr_float& S       = spdens(spid.metal(MetalSpLUT::S), 0);
  const gr_float& Fe      = spdens(spid.metal(MetalSpLUT::Fe), 0);
  const gr_float& MgSiO3  = spdens(spid.dust(DustSpLUT::MgSiO3_dust), 0);
  const gr_float& AC      = spdens(spid.dust(DustSpLUT::AC_dust), 0);
  const gr_float& SiM     = spdens(spid.dust(DustSpLUT::SiM_dust), 0);
  const gr_float& FeM     = spdens(spid.dust(DustSpLUT::FeM_dust), 0);
  const gr_float& Mg2SiO4 = spdens(spid.dust(DustSpLUT::Mg2SiO4_dust), 0);
  const gr_float& Fe3O4   = spdens(spid.dust(DustSpLUT::Fe3O4_dust), 0);
  const gr_float& SiO2D   = spdens(spid.dust(DustSpLUT::SiO2_dust), 0);
  const gr_float& MgO     = spdens(spid.dust(DustSpLUT::MgO_dust), 0);
  const gr_float& FeS     = spdens(spid.dust(DustSpLUT::FeS_dust), 0);
  const gr_float& Al2O3   = spdens(spid.dust(DustSpLUT::Al2O3_dust), 0);
  const gr_float& reforg  = spdens(spid.dust(DustSpLUT::ref_org_dust), 0);
  const gr_float& volorg  = spdens(spid.dust(DustSpLUT::vol_org_dust), 0);
  const gr_float& H2Oice  = spdens(spid.dust(DustSpLUT::H2O_ice_dust), 0);

  const double* const* kcol_buf = FullRxnRateBuf_kcol_bufs(&rxn_rate_buf);
  const double* const* kph_buf = FullRxnRateBuf_kph_bufs(&rxn_rate_buf);
  const double* const* grain_growth_rates =
    FullRxnRateBuf_grain_growth_bufs(&rxn_rate_buf);
  const double* h2dust = FullRxnRateBuf_h2dust(&rxn_rate_buf);
  // flag for if Gen Chiaki's dust model is enabled
  // -> to simplify bookkeeping in other parts of the code, we ignore whether
  //    the grain_growth flag was enabled, (the logic for setting the relevant
  //    reaction rates, explicitly sets the rates to 0 if disabled)
  const bool chiaki_model_dust_evolution = my_chemistry->dust_chemistry == 2;

  double scoef, acoef;

  // A) the 6-species integrator
  if (my_chemistry->primordial_chemistry == 1)  {




    // 1) HI

    scoef  = kcol_buf[CollisionalRxnLUT::k2][0]   *HII       *de;
    acoef  = kcol_buf[CollisionalRxnLUT::k1][0]   *de
           + kcol_buf[CollisionalRxnLUT::k57][0]   *HI
           + kcol_buf[CollisionalRxnLUT::k58][0]   *HeI       /4.
           + kph_buf[PhotoRxnLUT::k24][0];
    deriv(spid.prim(PrimordialSpLUT::HI), 0) =
        deriv(spid.prim(PrimordialSpLUT::HI), 0) +
        (scoef - acoef * HI);














    // 2) HII
    scoef  = kcol_buf[CollisionalRxnLUT::k1][0]   *HI    *de
           + kcol_buf[CollisionalRxnLUT::k57][0]   *HI    *HI
           + kcol_buf[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
           + kph_buf[PhotoRxnLUT::k24][0]   *HI;
    acoef  = kcol_buf[CollisionalRxnLUT::k2][0]   *de;
    deriv(spid.prim(PrimordialSpLUT::HII), 0) =
        deriv(spid.prim(PrimordialSpLUT::HII), 0) +
        (scoef - acoef * HII);
















    // 3) Electron density

    scoef = 0.
               + kcol_buf[CollisionalRxnLUT::k57][0]   *HI    *HI
               + kcol_buf[CollisionalRxnLUT::k58][0]   *HI    *HeI       /4.
               + kph_buf[PhotoRxnLUT::k24][0]   *HI
               + kph_buf[PhotoRxnLUT::k25][0]   *HeII       /4.
               + kph_buf[PhotoRxnLUT::k26][0]   *HeI       /4.;



    acoef = -(kcol_buf[CollisionalRxnLUT::k1][0]   *HI             - kcol_buf[CollisionalRxnLUT::k2][0]   *HII
            + kcol_buf[CollisionalRxnLUT::k3][0]   *HeI       /4. -
         kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            + kcol_buf[CollisionalRxnLUT::k5][0]   *HeII       /4. -
         kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       /4.);
    deriv(spid.prim(PrimordialSpLUT::e), 0) =
        deriv(spid.prim(PrimordialSpLUT::e), 0) +
        (scoef - acoef * de);





  }

  // --- (B) Do helium chemistry in any case: (for all ispecies values) ---




  // 4) HeI

  scoef  = kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       *de;
  acoef  = kcol_buf[CollisionalRxnLUT::k3][0]   *de
               + kph_buf[PhotoRxnLUT::k26][0];

  if (my_chemistry->primordial_chemistry > 3)  {
    scoef = scoef +  4. * ( 0.
        + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
        + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.
       );
    acoef = acoef
        + kcol_buf[CollisionalRxnLUT::k148][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k149][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k150][0]    *  H2II        /  2.;
  }
  deriv(spid.prim(PrimordialSpLUT::HeI), 0) =
      deriv(spid.prim(PrimordialSpLUT::HeI), 0) +
      (scoef - acoef * HeI);


  // 5) HeII

  scoef  = kcol_buf[CollisionalRxnLUT::k3][0]   *HeI    *de
         + kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       *de
         + kph_buf[PhotoRxnLUT::k26][0]   *HeI;

  acoef  = kcol_buf[CollisionalRxnLUT::k4][0]   *de        + kcol_buf[CollisionalRxnLUT::k5][0]   *de
         + kph_buf[PhotoRxnLUT::k25][0];

  if (my_chemistry->primordial_chemistry > 3)  {
    acoef = acoef
        + kcol_buf[CollisionalRxnLUT::k151][0]    *    HI;
  }
  deriv(spid.prim(PrimordialSpLUT::HeII), 0) =
      deriv(spid.prim(PrimordialSpLUT::HeII), 0) +
      (scoef - acoef * HeII);


  // 6) HeIII

  scoef   = kcol_buf[CollisionalRxnLUT::k5][0]   *HeII    *de
          + kph_buf[PhotoRxnLUT::k25][0]   *HeII;
  acoef   = kcol_buf[CollisionalRxnLUT::k6][0]   *de;
  deriv(spid.prim(PrimordialSpLUT::HeIII), 0) =
      deriv(spid.prim(PrimordialSpLUT::HeIII), 0) +
      (scoef - acoef * HeIII);





  // --- (C) Now do extra 3-species for molecular hydrogen ---

  if (my_chemistry->primordial_chemistry > 1)  {

    // First, do HI/HII with molecular hydrogen terms




    // 1) HI
    scoef  =      kcol_buf[CollisionalRxnLUT::k2][0]    * HII        * de
           + 2.*kcol_buf[CollisionalRxnLUT::k13][0]   * HI         * H2I       /2.
           +      kcol_buf[CollisionalRxnLUT::k11][0]   * HII        * H2I       /2.
           + 2.*kcol_buf[CollisionalRxnLUT::k12][0]   * de         * H2I       /2.
           +      kcol_buf[CollisionalRxnLUT::k14][0]   * HM         * de
           +      kcol_buf[CollisionalRxnLUT::k15][0]   * HM         * HI
           + 2.*kcol_buf[CollisionalRxnLUT::k16][0]   * HM         * HII
           + 2.*kcol_buf[CollisionalRxnLUT::k18][0]   * H2II       * de       /2.
           +      kcol_buf[CollisionalRxnLUT::k19][0]   * H2II       * HM       /2.
           + 2.*kph_buf[PhotoRxnLUT::k31][0]      * H2I       /2.;

    acoef  =      kcol_buf[CollisionalRxnLUT::k1][0]    * de
           +      kcol_buf[CollisionalRxnLUT::k7][0]    * de
           +      kcol_buf[CollisionalRxnLUT::k8][0]    * HM
           +      kcol_buf[CollisionalRxnLUT::k9][0]    * HII
           +      kcol_buf[CollisionalRxnLUT::k10][0]   * H2II       /2.
           + 2.*kcol_buf[CollisionalRxnLUT::k22][0]   * std::pow(HI       ,2)
           +      kcol_buf[CollisionalRxnLUT::k57][0]   * HI
           +      kcol_buf[CollisionalRxnLUT::k58][0]   * HeI       /4.
           + kph_buf[PhotoRxnLUT::k24][0];

    if (my_chemistry->primordial_chemistry > 2)  {
        scoef = scoef
          + kph_buf[PhotoRxnLUT::kdissHDI][0] * HDI       /3.0;
    }
    if ((my_chemistry->metal_chemistry == 1)  &&
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
        + kph_buf[PhotoRxnLUT::kdissOH][0]    * OHI        /17.0
        + kph_buf[PhotoRxnLUT::kdissH2O][0]   *H2OI       /18.0;
    }

    if (anydust != MASK_FALSE)  {
      if(itmask_metal[0] != MASK_FALSE   )  {
        acoef = acoef + 2. * h2dust[0]    * rhoH[0];
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k50][0]    * HII        * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k54][0]    * H2I        * DI         / 4.;
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k51][0]    * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k55][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM        *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
          + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        *    de        /  5.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        /  3.
          + kcol_buf[CollisionalRxnLUT::k130][0]    *   DII        /  2.
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k151][0]    *  HeII        /  4.
          + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        *   CHI       / 26.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        *   OHI        / 34.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *   OHI        *   COI        / 476.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        *   OHI        / 204.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        *   CHI       / 208.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        *   OHI        / 272.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII        *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII        *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII        *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        *   OHI        / 204.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 9.5
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        *   OHI        / 476.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        *   OHI        / 748.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *   CHI       / 13.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *  CH2I       / 14.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        / 28.;
    }
    deriv(spid.prim(PrimordialSpLUT::HI), 0) =
        deriv(spid.prim(PrimordialSpLUT::HI), 0) +
        (scoef - acoef * HI);













    // 2) HII

    scoef  =    kcol_buf[CollisionalRxnLUT::k1][0]     * HI        * de
           +    kcol_buf[CollisionalRxnLUT::k10][0]    * H2II       *HI       /2.
           +    kcol_buf[CollisionalRxnLUT::k57][0]    * HI        * HI
           +    kcol_buf[CollisionalRxnLUT::k58][0]    * HI        * HeI       /4.
           + kph_buf[PhotoRxnLUT::k24][0]   *HI;

    acoef  =    kcol_buf[CollisionalRxnLUT::k2][0]     * de
           +    kcol_buf[CollisionalRxnLUT::k9][0]     * HI
           +    kcol_buf[CollisionalRxnLUT::k11][0]    * H2I       /2.
           +    kcol_buf[CollisionalRxnLUT::k16][0]    * HM
           +    kcol_buf[CollisionalRxnLUT::k17][0]    * HM;
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k51][0]    * HI         * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k52][0]    * H2I        * DII        / 4.;
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k50][0]    * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k53][0]    * HDI        / 3.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k129][0]    *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k148][0]    *   HeI        /  4.
          + kcol_buf[CollisionalRxnLUT::k149][0]    *   HeI        /  4.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   O2I        / 32.;
    }
    deriv(spid.prim(PrimordialSpLUT::HII), 0) =
        deriv(spid.prim(PrimordialSpLUT::HII), 0) +
        (scoef - acoef * HII);

    
    // 3) electrons:

    scoef =   kcol_buf[CollisionalRxnLUT::k8][0]    * HM        * HI
           +  kcol_buf[CollisionalRxnLUT::k15][0]   * HM        * HI
           +  kcol_buf[CollisionalRxnLUT::k17][0]   * HM        * HII
           +  kcol_buf[CollisionalRxnLUT::k57][0]   * HI        * HI
           +  kcol_buf[CollisionalRxnLUT::k58][0]   * HI        * HeI       /4.
    // 
           + kph_buf[PhotoRxnLUT::k24][0]   *HI
           + kph_buf[PhotoRxnLUT::k25][0]   *HeII    /4.
           + kph_buf[PhotoRxnLUT::k26][0]   *HeI    /4.;

    if ((my_chemistry->metal_chemistry == 1)  &&
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef
        + kph_buf[PhotoRxnLUT::kphCI][0]   * CI       /12.0
        + kph_buf[PhotoRxnLUT::kphOI][0]   * OI       /16.0;
    }

    acoef = - (kcol_buf[CollisionalRxnLUT::k1][0]    *HI           - kcol_buf[CollisionalRxnLUT::k2][0]   *HII
            +  kcol_buf[CollisionalRxnLUT::k3][0]    *HeI       /4. -
         kcol_buf[CollisionalRxnLUT::k6][0]   *HeIII       /4.
            +  kcol_buf[CollisionalRxnLUT::k5][0]    *HeII       /4. -
         kcol_buf[CollisionalRxnLUT::k4][0]   *HeII       /4.
            +  kcol_buf[CollisionalRxnLUT::k14][0]   *HM
            -  kcol_buf[CollisionalRxnLUT::k7][0]    *HI
            -  kcol_buf[CollisionalRxnLUT::k18][0]   *H2II       /2.);
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef
            + kcol_buf[CollisionalRxnLUT::k56][0]    * DI         * HM        / 2.;
      acoef = acoef
            - kcol_buf[CollisionalRxnLUT::k1][0]     * DI         / 2.
            + kcol_buf[CollisionalRxnLUT::k2][0]     * DII        / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        /  3.
          + kcol_buf[CollisionalRxnLUT::k132][0]    *    DI        /  2.
          + kcol_buf[CollisionalRxnLUT::k153][0]    * HeHII        /  5.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        / 19.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        / 19.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *  O2II        / 32.;
    }
    deriv(spid.prim(PrimordialSpLUT::e), 0) =
        deriv(spid.prim(PrimordialSpLUT::e), 0) +
        (scoef - acoef * de);


    // 7) H2

    scoef = 2.*(kcol_buf[CollisionalRxnLUT::k8][0]     * HM          * HI
          +       kcol_buf[CollisionalRxnLUT::k10][0]    * H2II        * HI       /2.
          +       kcol_buf[CollisionalRxnLUT::k19][0]    * H2II        * HM       /2.
          +       kcol_buf[CollisionalRxnLUT::k22][0]    * HI        * std::pow((HI       ),2.));
    acoef = ( kcol_buf[CollisionalRxnLUT::k13][0]   *HI        + kcol_buf[CollisionalRxnLUT::k11][0]   *HII
            + kcol_buf[CollisionalRxnLUT::k12][0]   *de        )
            + kph_buf[PhotoRxnLUT::k29][0]    + kph_buf[PhotoRxnLUT::k31][0];

    if (anydust != MASK_FALSE)  {
      if(itmask_metal[0] != MASK_FALSE   )  {
        scoef = scoef + 2. * h2dust[0]    *
             HI        * rhoH[0];
      }
    }
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      scoef = scoef + 2. * (
              kcol_buf[CollisionalRxnLUT::k53][0]    * HDI        * HII        / 3.
            + kcol_buf[CollisionalRxnLUT::k55][0]    * HDI        * HI         / 3.
               );
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k52][0]    * DII        / 2.
            + kcol_buf[CollisionalRxnLUT::k54][0]    * DI         / 2.;
    }

    if ((my_chemistry->metal_chemistry == 1)  && 
        (itmask_metal[0] != MASK_FALSE))  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI        *   CHI       / 13.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI        *  CH2I       / 14.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI        *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI        *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
         );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   CHI       / 13.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        / 16.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        / 17.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        / 18.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *    CI        / 12.;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 0)  {
          scoef = scoef + 2. *
                grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      * 2.;

        }
        if (my_chemistry->dust_species > 1)  {
          scoef = scoef + 2. * (
                grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 3.
              + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       * 4.
              + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]
              + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       * 3.
            );
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      / H2I        * 2. * 2.;
        }
      }
    }
    deriv(spid.prim(PrimordialSpLUT::H2I), 0) =
        deriv(spid.prim(PrimordialSpLUT::H2I), 0) +
        (scoef - acoef * H2I);


    // 8) H-

    scoef = kcol_buf[CollisionalRxnLUT::k7][0]    * HI        * de;
    acoef = (kcol_buf[CollisionalRxnLUT::k8][0]     + kcol_buf[CollisionalRxnLUT::k15][0]   )  * HI        +
            (kcol_buf[CollisionalRxnLUT::k16][0]    + kcol_buf[CollisionalRxnLUT::k17][0]   )  * HII        +
            kcol_buf[CollisionalRxnLUT::k14][0]    * de        + kcol_buf[CollisionalRxnLUT::k19][0]    * H2II       /2.0f +
            kph_buf[PhotoRxnLUT::k27][0];
    // contribution of minor species
    if (my_chemistry->primordial_chemistry > 2)  {
      acoef = acoef
            + kcol_buf[CollisionalRxnLUT::k56][0]    * DI         / 2.;
    }

    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.;
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    DI        /  2.;
    }
    deriv(spid.prim(PrimordialSpLUT::HM), 0) =
        deriv(spid.prim(PrimordialSpLUT::HM), 0) +
        (scoef - acoef * HM);



    // 9) H2+

    scoef =    2.*( kcol_buf[CollisionalRxnLUT::k9][0]    *HI    *HII
                  +   kcol_buf[CollisionalRxnLUT::k11][0]   *H2I    /2.*HII
                  +   kcol_buf[CollisionalRxnLUT::k17][0]   *HM    *HII
                  + kph_buf[PhotoRxnLUT::k29][0]   *H2I    /2.
                  );
    acoef =         kcol_buf[CollisionalRxnLUT::k10][0]   *HI     + kcol_buf[CollisionalRxnLUT::k18][0]   *de
                  + kcol_buf[CollisionalRxnLUT::k19][0]   *HM
                  + (kph_buf[PhotoRxnLUT::k28][0]   +kph_buf[PhotoRxnLUT::k30][0]   );
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k152][0]    * HeHII        *    HI        /  5.
         );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        /  4.;
    }
    deriv(spid.prim(PrimordialSpLUT::H2II), 0) =
        deriv(spid.prim(PrimordialSpLUT::H2II), 0) +
        (scoef - acoef * H2II);









  }

  // --- (D) Now do extra 3-species for molecular HD ---
  if (my_chemistry->primordial_chemistry > 2)  {


    
    // 1) DI
    scoef =   (       kcol_buf[CollisionalRxnLUT::k2][0]    * DII        * de
               +      kcol_buf[CollisionalRxnLUT::k51][0]   * DII        * HI
               + 2.*kcol_buf[CollisionalRxnLUT::k55][0]   * HDI        *
            HI       /3.
               )
               + 2. * kph_buf[PhotoRxnLUT::kdissHDI][0]  *   HDI /3.;
    acoef  =    kcol_buf[CollisionalRxnLUT::k1][0]    * de
           +    kcol_buf[CollisionalRxnLUT::k50][0]    * HII
           +    kcol_buf[CollisionalRxnLUT::k54][0]    * H2I       /2.
           +    kcol_buf[CollisionalRxnLUT::k56][0]    * HM
           + kph_buf[PhotoRxnLUT::k24][0];
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  2. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k131][0]    *  HDII        *    de        /  3.
          + kcol_buf[CollisionalRxnLUT::k133][0]    *   DII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII        *    DM        /  2.
          + kcol_buf[CollisionalRxnLUT::k136][0]    *    DM        *    HI        /  2.
          );
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k129][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::k132][0]    *    de
          + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM;
    }
    deriv(spid.prim(PrimordialSpLUT::DI), 0) =
        deriv(spid.prim(PrimordialSpLUT::DI), 0) +
        (scoef - acoef * DI);
                                                    

    // 2) DII
    scoef =   (   kcol_buf[CollisionalRxnLUT::k1][0]     * DI        * de
          +       kcol_buf[CollisionalRxnLUT::k50][0]    * HII       * DI
          +  2.*kcol_buf[CollisionalRxnLUT::k53][0]    * HII       * HDI       /3.
          )
          + kph_buf[PhotoRxnLUT::k24][0]   *DI;
    acoef = 0.;
    // ! initialize GC202002
    acoef =    kcol_buf[CollisionalRxnLUT::k2][0]     * de
          +    kcol_buf[CollisionalRxnLUT::k51][0]    * HI
          +    kcol_buf[CollisionalRxnLUT::k52][0]    * H2I       /2.;
    if (my_chemistry->primordial_chemistry > 3)  {
      acoef = acoef
          + kcol_buf[CollisionalRxnLUT::k130][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::k133][0]    *    DM        /  2.;
    }
    deriv(spid.prim(PrimordialSpLUT::DII), 0) =
        deriv(spid.prim(PrimordialSpLUT::DII), 0) +
        (scoef - acoef * DII);


    // 3) HDI
    scoef = 3.*(kcol_buf[CollisionalRxnLUT::k52][0]    * DII       *
         H2I       /2./2.
         + kcol_buf[CollisionalRxnLUT::k54][0]    * DI        * H2I       /2./2.
    // !   &           + 2._DKIND*kcol_buf[CollisionalRxnLUT::k56][0]    * DI        * HM       /2._DKIND
    //- ! corrected by GC202005
         +          kcol_buf[CollisionalRxnLUT::k56][0]    * DI        * HM       /2.
               );
    acoef  =    kcol_buf[CollisionalRxnLUT::k53][0]    * HII
           +    kcol_buf[CollisionalRxnLUT::k55][0]    * HI
           +  kph_buf[PhotoRxnLUT::kdissHDI][0];
    if (my_chemistry->primordial_chemistry > 3)  {
      scoef = scoef +  3. * ( 0.
          + kcol_buf[CollisionalRxnLUT::k125][0]    *  HDII        *    HI        /  3.
          + kcol_buf[CollisionalRxnLUT::k137][0]    *    DM        *    HI        /  2.
          );
    }
    deriv(spid.prim(PrimordialSpLUT::HDI), 0) =
        deriv(spid.prim(PrimordialSpLUT::HDI), 0) +
        (scoef - acoef * HDI);




  }

  // --- (D2) Now do extra 3-species for minor primordial species ---
  if (my_chemistry->primordial_chemistry > 3)  {



    // 1) DM

    scoef =
          kcol_buf[CollisionalRxnLUT::k132][0]    *    DI        *    de
        + kcol_buf[CollisionalRxnLUT::k135][0]    *    HM        *    DI;
    acoef =
          kcol_buf[CollisionalRxnLUT::k133][0]    *   DII        /  2.
        + kcol_buf[CollisionalRxnLUT::k134][0]    *   HII
        + kcol_buf[CollisionalRxnLUT::k136][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k137][0]    *    HI;

    deriv(spid.prim(PrimordialSpLUT::DM), 0) =
        deriv(spid.prim(PrimordialSpLUT::DM), 0) +
        (scoef - acoef * DM);


    // 2) HDII

    scoef = 3. * (
          kcol_buf[CollisionalRxnLUT::k129][0]    *    DI        *   HII        /  2.
        + kcol_buf[CollisionalRxnLUT::k130][0]    *   DII        *    HI        /  2.
       );
    acoef =
          kcol_buf[CollisionalRxnLUT::k125][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k131][0]    *    de;

    deriv(spid.prim(PrimordialSpLUT::HDII), 0) =
        deriv(spid.prim(PrimordialSpLUT::HDII), 0) +
        (scoef - acoef * HDII);


    // 3) HeHII

    scoef = 5. * (
          kcol_buf[CollisionalRxnLUT::k148][0]    *   HeI        *   HII        /  4.
        + kcol_buf[CollisionalRxnLUT::k149][0]    *   HeI        *   HII        /  4.
        + kcol_buf[CollisionalRxnLUT::k150][0]    *   HeI        *  H2II        /  8.
        + kcol_buf[CollisionalRxnLUT::k151][0]    *  HeII        *    HI        /  4.
       );
    acoef =
          kcol_buf[CollisionalRxnLUT::k152][0]    *    HI
        + kcol_buf[CollisionalRxnLUT::k153][0]    *    de;

    deriv(spid.prim(PrimordialSpLUT::HeHII), 0) =
        deriv(spid.prim(PrimordialSpLUT::HeHII), 0) +
        (scoef - acoef * HeHII);




  }

  // --- (D3) Now do metal species ---
  if (my_chemistry->metal_chemistry == 1)  {

    if (itmask_metal[0] != MASK_FALSE   )  {

      // ***** CI **********
      scoef = 0. + 12. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI        *   CHI       / 13.
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *   CII        *    de        / 12.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *   H2I        /  2.
          + kph_buf[PhotoRxnLUT::kphCI][0];
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::AC_dust][0]          / CI        * 12.;
        }
      }
      scoef = scoef + 12. *
          kph_buf[PhotoRxnLUT::kdissCO][0]    *COI        /28.0;

      deriv(spid.metal(MetalSpLUT::CI), 0) =
          deriv(spid.metal(MetalSpLUT::CI), 0) +
          (scoef - acoef * CI);



      // ***** CII **********
      scoef = 0. + 12. * ( 0.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz44][0]    *    de;
      scoef = scoef
          + kph_buf[PhotoRxnLUT::kphCI][0]    * CI;

      deriv(spid.metal(MetalSpLUT::CII), 0) =
          deriv(spid.metal(MetalSpLUT::CII), 0) +
          (scoef - acoef * CII);



      // ***** COI **********
      scoef = 0. + 28. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        *   OHI        / 204.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        *   O2I        / 384.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        *   CHI       / 208.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        *   O2I        / 384.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *  COII        *    HI        / 28.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *   OHI        / 17.;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      /COI        * 17. * 0.5
          + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      /COI        * 17.;
        }
      }
      acoef = acoef
          + kph_buf[PhotoRxnLUT::kdissCO][0];

      deriv(spid.metal(MetalSpLUT::COI), 0) =
          deriv(spid.metal(MetalSpLUT::COI), 0) +
          (scoef - acoef * COI);



      // ***** CO2I *********
      scoef = 0. + 44. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *   OHI        *   COI        / 476.
         );
      acoef = 0.;

      deriv(spid.metal(MetalSpLUT::CO2I), 0) =
          deriv(spid.metal(MetalSpLUT::CO2I), 0) +
          (scoef - acoef * CO2I);



      // ***** OI **********
      scoef = 0. + 16. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI        *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI        *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *   OHI        *   OHI        / 289.
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        *   O2I        / 384.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *   OII        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *   OII        *    de        / 16.
          + kcol_buf[CollisionalRxnLUT::kz47][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *  O2II        *    de        / 16.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        *   O2I        / 896.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz31][0]    *    OI        / 8.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *   CHI       / 13.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *   OHI        / 17.
          + kph_buf[PhotoRxnLUT::kphOI][0];
      scoef = scoef + 16. *
          ( kph_buf[PhotoRxnLUT::kdissOH][0]    *OHI        /17.0
          + kph_buf[PhotoRxnLUT::kdissCO][0]    *COI        /28.0);

      deriv(spid.metal(MetalSpLUT::OI), 0) =
          deriv(spid.metal(MetalSpLUT::OI), 0) +
          (scoef - acoef * OI);



      // ***** OH *********
      scoef = 0. + 17. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI        *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI        *   O2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz21][0]    *    OI        *   H2I        / 32.
          + kcol_buf[CollisionalRxnLUT::kz30][0]    *    OI        *    HI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    * H2OII        *    de        / 18.
          + kcol_buf[CollisionalRxnLUT::kz49][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz17][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *   OHI        / 8.5
          + kcol_buf[CollisionalRxnLUT::kz26][0]    *   COI        / 28.
          + kcol_buf[CollisionalRxnLUT::kz28][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        / 28.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        / 44.;
      acoef = acoef
          + kph_buf[PhotoRxnLUT::kdissOH][0];
      scoef = scoef + 17. *
          kph_buf[PhotoRxnLUT::kdissH2O][0]    *H2OI      /18.0;

      deriv(spid.metal(MetalSpLUT::OHI), 0) =
          deriv(spid.metal(MetalSpLUT::OHI), 0) +
          (scoef - acoef * OHI);



      // ***** H2O *********
      scoef = 0. + 18. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz24][0]    *   H2I        *   OHI        / 34.
          + kcol_buf[CollisionalRxnLUT::kz25][0]    *   OHI        *   OHI        / 289.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    * H3OII        *    de        / 19.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz18][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      /H2OI       * 18. * 2.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     /H2OI       * 18. * 3.
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       /H2OI       * 18. * 4.
          + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         /H2OI       * 18.
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       /H2OI       * 18. * 3.;
        }
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust][0]      /H2OI       * 18.;
        }
      }
      acoef = acoef
          + kph_buf[PhotoRxnLUT::kdissH2O][0];

      deriv(spid.metal(MetalSpLUT::H2OI), 0) =
          deriv(spid.metal(MetalSpLUT::H2OI), 0) +
          (scoef - acoef * H2OI);



      // ***** O2I *********
      scoef = 0. + 32. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz31][0]    *    OI        *    OI        / 256.
          + kcol_buf[CollisionalRxnLUT::kz33][0]    *    OI        *   OHI        / 272.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz19][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz29][0]    *    CI        / 12.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        / 12.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        / 28.;

      deriv(spid.metal(MetalSpLUT::O2I), 0) =
          deriv(spid.metal(MetalSpLUT::O2I), 0) +
          (scoef - acoef * O2I);



      // ***** SiI **********
      scoef = 0. + 28. * ( 0.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   O2I        / 32.;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::SiM_dust][0]         / SiI        * 28.;
        }
      }

      deriv(spid.metal(MetalSpLUT::SiI), 0) =
          deriv(spid.metal(MetalSpLUT::SiI), 0) +
          (scoef - acoef * SiI);



      // ***** SiOI **********
      scoef = 0. + 44. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz52][0]    *   SiI        *   OHI        / 476.
          + kcol_buf[CollisionalRxnLUT::kz53][0]    *   SiI        *   O2I        / 896.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *   OHI        / 17.;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 0)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      / SiOI        * 44.;
        }
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / SiOI        * 44.;
        }
      }

      deriv(spid.metal(MetalSpLUT::SiOI), 0) =
          deriv(spid.metal(MetalSpLUT::SiOI), 0) +
          (scoef - acoef * SiOI);



      // ***** SiO2I **********
      scoef = 0. + 60. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz54][0]    *  SiOI        *   OHI        / 748.
         );
      acoef = 0.;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 1)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust][0]       / SiO2I        * 60.;
        }
      }

      deriv(spid.metal(MetalSpLUT::SiO2I), 0) =
          deriv(spid.metal(MetalSpLUT::SiO2I), 0) +
          (scoef - acoef * SiO2I);



      // ***** CHI *********
      scoef = 0. + 13. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI        *  CH2I       / 14.
          + kcol_buf[CollisionalRxnLUT::kz20][0]    *    CI        *   H2I        / 24.
          + kcol_buf[CollisionalRxnLUT::kz27][0]    *    CI        *    HI        / 12.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz15][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz32][0]    *    OI        / 16.;

      deriv(spid.metal(MetalSpLUT::CHI), 0) =
          deriv(spid.metal(MetalSpLUT::CHI), 0) +
          (scoef - acoef * CHI);



      // ***** CH2I *********
      scoef = 0. + 14. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz23][0]    *   H2I        *   CHI       / 26.
          + kcol_buf[CollisionalRxnLUT::kz51][0]    *   H2I        *    CI        / 24.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz16][0]    *    HI;
      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 2)  {
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      /CH2I       * 14. * 0.5;
        }
      }

      deriv(spid.metal(MetalSpLUT::CH2I), 0) =
          deriv(spid.metal(MetalSpLUT::CH2I), 0) +
          (scoef - acoef * CH2I);



      // ***** COII **********
      scoef = 0. + 28. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz37][0]    *   CII        *   OHI        / 204.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz43][0]    *    HI;

      deriv(spid.metal(MetalSpLUT::COII), 0) =
          deriv(spid.metal(MetalSpLUT::COII), 0) +
          (scoef - acoef * COII);



      // ***** OII **********
      scoef = 0. + 16. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz22][0]    *   HII        *    OI        / 16.
          + kcol_buf[CollisionalRxnLUT::kz38][0]    *   CII        *   O2I        / 384.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz39][0]    *    HI
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz45][0]    *    de;
      scoef = scoef
          + kph_buf[PhotoRxnLUT::kphOI][0]    * OI;

      deriv(spid.metal(MetalSpLUT::OII), 0) =
          deriv(spid.metal(MetalSpLUT::OII), 0) +
          (scoef - acoef * OII);



      // ***** OHII **********
      scoef = 0. + 17. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz34][0]    *   HII        *   OHI        / 17.
          + kcol_buf[CollisionalRxnLUT::kz40][0]    *   OII        *   H2I        / 32.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *   H2I        /  2.;

      deriv(spid.metal(MetalSpLUT::OHII), 0) =
          deriv(spid.metal(MetalSpLUT::OHII), 0) +
          (scoef - acoef * OHII);



      // ***** H2OII **********
      scoef = 0. + 18. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz35][0]    *   HII        *  H2OI       / 18.
          + kcol_buf[CollisionalRxnLUT::kz41][0]    *  OHII        *   H2I        / 34.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    *   H2I        /  2.
          + kcol_buf[CollisionalRxnLUT::kz46][0]    *    de
          + kcol_buf[CollisionalRxnLUT::kz47][0]    *    de;

      deriv(spid.metal(MetalSpLUT::H2OII), 0) =
          deriv(spid.metal(MetalSpLUT::H2OII), 0) +
          (scoef - acoef * H2OII);



      // ***** H3OII **********
      scoef = 0. + 19. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz42][0]    * H2OII        *   H2I        / 36.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz48][0]    *    de
          + kcol_buf[CollisionalRxnLUT::kz49][0]    *    de;

      deriv(spid.metal(MetalSpLUT::H3OII), 0) =
          deriv(spid.metal(MetalSpLUT::H3OII), 0) +
          (scoef - acoef * H3OII);



      // ***** O2II **********
      scoef = 0. + 32. * ( 0.
          + kcol_buf[CollisionalRxnLUT::kz36][0]    *   HII        *   O2I        / 32.
         );
      acoef = 0.
          + kcol_buf[CollisionalRxnLUT::kz50][0]    *    de;

      deriv(spid.metal(MetalSpLUT::O2II), 0) =
          deriv(spid.metal(MetalSpLUT::O2II), 0) +
          (scoef - acoef * O2II);



      if (chiaki_model_dust_evolution)  {
        if (my_chemistry->dust_species > 0)  {
          // ***** Mg **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      / Mg        * 24.;
          if (my_chemistry->dust_species > 1)  {
            acoef = acoef
            + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     / Mg        * 24. * 2.
            + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         / Mg        * 24.;
          }

          deriv(spid.metal(MetalSpLUT::Mg), 0) =
              deriv(spid.metal(MetalSpLUT::Mg), 0) +
              (scoef - acoef * Mg);


        }

        if (my_chemistry->dust_species > 1)  {
          // ***** Al **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       / Al        * 27. * 2.;

          deriv(spid.metal(MetalSpLUT::Al), 0) =
              deriv(spid.metal(MetalSpLUT::Al), 0) +
              (scoef - acoef * Al);



          // ***** S  **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         / S        * 32.;

          deriv(spid.metal(MetalSpLUT::S), 0) =
              deriv(spid.metal(MetalSpLUT::S), 0) +
              (scoef - acoef * S);



          // ***** Fe **********
          scoef = 0.;
          acoef = 0.;
          acoef = acoef
          + grain_growth_rates[OnlyGrainSpLUT::FeM_dust][0]         / Fe        * 56.
          + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       / Fe        * 56. * 3.
          + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         / Fe        * 56.;

          deriv(spid.metal(MetalSpLUT::Fe), 0) =
              deriv(spid.metal(MetalSpLUT::Fe), 0) +
              (scoef - acoef * Fe);


        }
      }

    }

  }

  // --- (D4) Now do dust species ---
  if (chiaki_model_dust_evolution)  {

    if (itmask_metal[0] != MASK_FALSE   )  {

      if (my_chemistry->dust_species > 0)  {
        // ***** MgSiO3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][0]      * 100.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::MgSiO3_dust), 0) =
            deriv(spid.dust(DustSpLUT::MgSiO3_dust), 0) +
            (scoef - acoef * MgSiO3);



        // ***** AC **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::AC_dust][0]          * 12.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::AC_dust), 0) =
            deriv(spid.dust(DustSpLUT::AC_dust), 0) +
            (scoef - acoef * AC);


      }

      if (my_chemistry->dust_species > 1)  {
        // ***** SiM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::SiM_dust][0]         * 28.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::SiM_dust), 0) =
            deriv(spid.dust(DustSpLUT::SiM_dust), 0) +
            (scoef - acoef * SiM);



        // ***** FeM **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::FeM_dust][0]         * 56.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::FeM_dust), 0) =
            deriv(spid.dust(DustSpLUT::FeM_dust), 0) +
            (scoef - acoef * FeM);



        // ***** Mg2SiO4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Mg2SiO4_dust][0]     * 140.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::Mg2SiO4_dust), 0) =
            deriv(spid.dust(DustSpLUT::Mg2SiO4_dust), 0) +
            (scoef - acoef * Mg2SiO4);



        // ***** Fe3O4 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Fe3O4_dust][0]       * 232.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::Fe3O4_dust), 0) =
            deriv(spid.dust(DustSpLUT::Fe3O4_dust), 0) +
            (scoef - acoef * Fe3O4);



        // ***** SiO2D **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::SiO2_dust][0]       * 60.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::SiO2_dust), 0) =
            deriv(spid.dust(DustSpLUT::SiO2_dust), 0) +
            (scoef - acoef * SiO2D);



        // ***** MgO **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::MgO_dust][0]         * 40.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::MgO_dust), 0) =
            deriv(spid.dust(DustSpLUT::MgO_dust), 0) +
            (scoef - acoef * MgO);



        // ***** FeS **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::FeS_dust][0]         * 88.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::FeS_dust), 0) =
            deriv(spid.dust(DustSpLUT::FeS_dust), 0) +
            (scoef - acoef * FeS);



        // ***** Al2O3 **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::Al2O3_dust][0]       * 102.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::Al2O3_dust), 0) =
            deriv(spid.dust(DustSpLUT::Al2O3_dust), 0) +
            (scoef - acoef * Al2O3);


      }

      if (my_chemistry->dust_species > 2)  {
        // ***** reforg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::ref_org_dust][0]      * 22.68;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::ref_org_dust), 0) =
            deriv(spid.dust(DustSpLUT::ref_org_dust), 0) +
            (scoef - acoef * reforg);



        // ***** volorg **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::vol_org_dust][0]      * 32.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::vol_org_dust), 0) =
            deriv(spid.dust(DustSpLUT::vol_org_dust), 0) +
            (scoef - acoef * volorg);



        // ***** H2Oice **********
        scoef = 0.;
        scoef = scoef
        + grain_growth_rates[OnlyGrainSpLUT::H2O_ice_dust][0]      * 18.;
        acoef = 0.;

        deriv(spid.dust(DustSpLUT::H2O_ice_dust), 0) =
            deriv(spid.dust(DustSpLUT::H2O_ice_dust), 0) +
            (scoef - acoef * H2Oice);


      }

    }

  }

}




} // namespace grackle::impl::chemistry

#endif /* CHEMISTRY_SOLVER_FUNCS_HPP */
