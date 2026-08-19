//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the @ref DustSolver type
///
//===----------------------------------------------------------------------===//

#include <cfloat>  // DBL_MAX

#include "grackle.h"
#include "dust/grain_species_info.hpp"
#include "dust/multi_grain_species/calc_grain_size_increment_1d.hpp"
#include "dust/solver.hpp"
#include "interpolate.hpp"
#include "opaque_storage.hpp"
#include "utils-cpp.hpp"
#include "utils-field.hpp"

namespace GRIMPL_NAMESPACE_DECL {

// Some ideas for refactoring
// - we probably want to break the following function up into its constituent
//   parts.
// - If nothing else, it would be nice to more cleanly separate logic for the
//   different dust models (i.e. variants of the classic 1-field model and the
//   multi-grain-species model)
// - we may also want to give some thought to possibly grouping subsets of the
//   arguments that are only used for certain dust models.

void DustSolver::lookup_dust_rates1d(
    IndexRange idx_range, const double* tdust, const double* dust2gas,
    double dom, const gr_mask_type* itmask_metal, chemistry_data* my_chemistry,
    chemistry_data_storage* my_rates, grackle_field_data* my_fields,
    GrainSpeciesCollection grain_temperatures,
    LnTLinInterpBuf logTlininterp_buf, FullRxnRateBuf rxn_rate_buf,
    grackle::impl::InternalDustPropBuf internal_dust_prop_scratch_buf) {
  // TODO: get rid of dlogtem argument!

  const double dlogTdust =
      (std::log(my_chemistry->DustTemperatureEnd) -
       std::log(my_chemistry->DustTemperatureStart)) /
      (double)(my_chemistry->NumberOfDustTemperatureBins - 1);

  // we should probably enforce the following at initialization!
  GRIMPL_REQUIRE(my_chemistry->dust_species >= 0, "sanity-check!");

  // h2dust is the buffer that gets filled with the rate for forming
  // molecular hydrogen on dust grains
  double* h2dust = FullRxnRateBuf_h2dust(&rxn_rate_buf);

  if (my_chemistry->dust_species == 0) {
    // in this branch, we are just tracking a single generic dust field

    const double logTdust_start = std::log(my_chemistry->DustTemperatureStart);
    const double logTdust_end = std::log(my_chemistry->DustTemperatureEnd);
    // construct a view the h2dust interpolation table
    FortranView<double**> h2dusta(my_rates->h2dust,
                                  my_chemistry->NumberOfTemperatureBins,
                                  my_chemistry->NumberOfDustTemperatureBins);

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        // Assume dust melts at Tdust > DustTemperatureEnd, in the context of
        // computing the H2 formation rate
        //
        // important: at the time of writing, when using a generic dust
        // density field (i.e. my_chemistry->use_dust_density_field == 1),
        // we don't mutate that density field. This contrasts with Grackle's
        // behavior when tracking dust species fields.
        if (tdust[i] > my_chemistry->DustTemperatureEnd) {
          h2dust[i] = tiny8;
        } else {
          // Get log dust temperature

          double logTdust = grackle::impl::clamp(std::log(tdust[i]),
                                                 logTdust_start, logTdust_end);

          // Find index into table and precompute interpolation values

          long long d_indixe = grackle::impl::clamp(
              (long long)((logTdust - logTdust_start) / dlogTdust) + 1LL, 1LL,
              (long long)(my_chemistry->NumberOfDustTemperatureBins) - 1LL);
          double d_t1 = (logTdust_start + (d_indixe - 1) * dlogTdust);
          double d_t2 = (logTdust_start + d_indixe * dlogTdust);
          double d_tdef = (logTdust - d_t1) / (d_t2 - d_t1);

          // Get rate from 2D interpolation

          double dusti1 =
              h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe - 1) +
              (h2dusta(logTlininterp_buf.indixe[i], d_indixe - 1) -
               h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe - 1)) *
                  logTlininterp_buf.tdef[i];
          double dusti2 = h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe) +
                          (h2dusta(logTlininterp_buf.indixe[i], d_indixe) -
                           h2dusta(logTlininterp_buf.indixe[i] - 1, d_indixe)) *
                              logTlininterp_buf.tdef[i];
          h2dust[i] = dusti1 + (dusti2 - dusti1) * d_tdef;

          // Multiply by dust to gas ratio

          h2dust[i] = h2dust[i] * dust2gas[i];
        }
      }
    }

  } else {  // my_chemistry->dust_species > 0
    const double dlogtem = common_1D_rate_table_lnT_step(*my_chemistry);

    // The Fortran version of this function implicitly assumed that the
    // following condition was satisfied when my_chemistry->dust_species > 0
    // (we are just making it more explicit)
    GRIMPL_REQUIRE(my_chemistry->metal_chemistry == 1, "sanity-check!");

    // Compute grain size increment
    calc_grain_size_increment_1d(dom, idx_range, itmask_metal, my_chemistry,
                                 my_rates->opaque_storage->grain_species_info,
                                 my_rates->opaque_storage->inject_pathway_props,
                                 my_fields, internal_dust_prop_scratch_buf);

    FortranView<const gr_float***> d(
        my_fields->density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    // the use of SpeciesLUTFieldAdaptor with dynamic indices is suboptimal
    // (it triggers indices). But, I think its ok here for 2 reasons:
    // 1. we aren't using it deep in a loop
    // 2. it makes the code significantly easier to manage! Plus, if somebody
    //    is using this configuration, speed is clearly not a huge priority
    //    (given all of the extra passive scalars that must be tracked)
    grackle::impl::SpeciesLUTFieldAdaptor field_data_adaptor{*my_fields};

    grackle::impl::GrainSpeciesInfo* gsp_info =
        my_rates->opaque_storage->grain_species_info;
    GRIMPL_REQUIRE(gsp_info != nullptr, "sanity check!");

    int n_grain_species = gsp_info->n_species;

    // Look-up rate for H2 formation on dust
    // -------------------------------------
    // in this branch, we are effectively computing:
    //   h2dust[i]
    //     = ∑ₛ coef_fnₛ(Tgas[i], Tdustₛ[i]) * grain_sigma_per_gas_massₛ[i]
    // where, the "s" subscript corresponds to a dust species.
    //
    // In practice, `coef_fnₛ` is implemented through 2D interpolation
    // of tabulated values

    // load properties of the interpolation table
    InterpGridProps& interp_props =
        my_rates->opaque_storage->h2dust_grain_interp_props;

    // load the tables that we are interpolating over
    //
    // TODO: remove these local variables
    // - strictly speaking, these local variables are **NOT** necessary even
    //   now. They just exist for expositionary purposes
    // - we should only delete these local variables **AFTER** the struct
    //   members holding the values have better names
    // - both of the following:
    //     1. the choice to **ONLY** use silicate & carbonaceous coefficients
    //        for computing H2 formation rates (when several species are
    //        neither silicates nor carbonaceous)
    //     2. the choice of which coefficients to user for each grain species
    //   all seem a little strange (perhaps even questionable)... We go into
    //   more detail in a comment when we initialize the grain species
    //   information
    const double* h2rate_silicate_coef_table = my_rates->h2dustS;
    const double* h2rate_carbonaceous_coef_table = my_rates->h2dustC;

    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      if (itmask_metal[i] != MASK_FALSE) {
        double summed_h2dust = 0.0;
        for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
          const double* coef_table =
              gsp_info->species_info[gsp_idx].h2dust_uses_carbonaceous_table
                  ? h2rate_carbonaceous_coef_table
                  : h2rate_silicate_coef_table;
          // take the natural log of the grain species's Temperature
          double logTdust = std::log(grain_temperatures.data[gsp_idx][i]);
          double coef = interpolate_2d(
              logTdust, logTlininterp_buf.logtem[i], interp_props.dimension,
              interp_props.parameters[0], dlogTdust, interp_props.parameters[1],
              dlogtem, interp_props.data_size, coef_table);
          double sigma_per_gas_mass =
              internal_dust_prop_scratch_buf.grain_sigma_per_gas_mass
                  .data[gsp_idx][i];
          summed_h2dust += coef * sigma_per_gas_mass;
        }
        h2dust[i] = summed_h2dust;
      }
    }

    // Compute net grain growth rates
    // ------------------------------
    double* const* grain_growth_rates =
        FullRxnRateBuf_grain_growth_bufs(&rxn_rate_buf);

    long long nratec_single_elem_arr[1] = {
        (long long)(my_chemistry->NumberOfTemperatureBins)};

    if (my_chemistry->grain_growth == 1) {
      // NOTE: an earlier version of this logic included a large block of
      //       commented code that appeared to calculate growth rates for
      //       the MgSiO3 dust grain species
      //  - the block included a comment stating that
      //    "Formulation from Nozawa et al. (2003, 2012)"
      //  - importantly, this logic totally deviated from the highly
      //    standard logic that we currently use for computing growth rates
      //    MgSiO3 and all other grain species

      // iterate over the grain species
      for (int gsp_idx = 0; gsp_idx < n_grain_species; gsp_idx++) {
        //  get the number of ingredients for the current grain species &
        //  the actual ingredient list
        int n_ingred = gsp_info->species_info[gsp_idx].n_growth_ingredients;
        const GrainGrowthIngredient* ingredient_l =
            gsp_info->species_info[gsp_idx].growth_ingredients;

        // preload views of each ingredient's current mass density &
        // precompute the divisor divided by the mass density
        FortranView<const gr_float***>
            ingred_view[grackle::impl::max_ingredients_per_grain_species];
        double ingred_divisor[grackle::impl::max_ingredients_per_grain_species];

        for (int ingred_idx = 0; ingred_idx < n_ingred; ingred_idx++) {
          const gr_float* ptr = field_data_adaptor.get_ptr_dynamic(
              ingredient_l[ingred_idx].species_idx);
          ingred_view[ingred_idx] = FortranView<const gr_float***>(
              ptr, my_fields->grid_dimension[0], my_fields->grid_dimension[1],
              my_fields->grid_dimension[2]);

          // load the particle mass for the current ingredient
          double mparticle_amu = ingredient_l[ingred_idx].mparticle_amu;
          // load the stoichiometric coefficient
          double coef = ingredient_l[ingred_idx].coef;
          ingred_divisor[ingred_idx] = std::pow(mparticle_amu, 1.5) * coef;
        }

        // preload the grain species's grain_sigma_per_gas_mass
        const double* grain_sigma_per_gas_mass =
            internal_dust_prop_scratch_buf.grain_sigma_per_gas_mass
                .data[gsp_idx];

        // compute the growth rate at each location in idx_range
        for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
          if (itmask_metal[i] != MASK_FALSE) {
            // calculate the factor that depends on the limiting reactant
            double limiting_factor = DBL_MAX;  // <- max finite value
            for (int ingred_idx = 0; ingred_idx < n_ingred; ingred_idx++) {
              limiting_factor = std::fmin(
                  limiting_factor,
                  ingred_view[ingred_idx](i, idx_range.j, idx_range.k) /
                      ingred_divisor[ingred_idx]);
            }
            // if ingredient_l is empty, set limiting_factor to 0
            limiting_factor *= (n_ingred > 0);

            // finally, we're ready to compute the rate
            double kd = interpolate_1d(
                logTlininterp_buf.logtem[i], nratec_single_elem_arr,
                interp_props.parameters[1], dlogtem, nratec_single_elem_arr[0],
                my_rates->grain_growth_rate);

            grain_growth_rates[OnlyGrainSpLUT::MgSiO3_dust][i] =
                kd * grain_sigma_per_gas_mass[i] *
                d(i, idx_range.j, idx_range.k) * limiting_factor;
          }
        }  // idx_range loop
      }  // n_grain_species loop
    }
  }
}

}  // namespace GRIMPL_NAMESPACE_DECL