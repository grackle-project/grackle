#ifndef DUST_GROWTH_AND_DESTRUCTION_HPP
#define DUST_GROWTH_AND_DESTRUCTION_HPP

#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "index_helper.h"
#include "internal_units.hpp"
#include "phys_constants.h"
#include "fortran_func_decls.h"

namespace grackle::impl {

// Calculates dust growth rates (accretion) onto grain surfaces.
// Stores the mass change dM for each cell in growth_dM array.
void dust_growth(chemistry_data* my_chemistry, grackle_field_data* my_fields,
                 InternalGrUnits internalu, IndexRange idx_range,
                 const gr_mask_type* itmask, const double* dt_value,
                 const double* t_gas,
                 double* growth_dM  // output: mass change rate for each cell
);

// Species-specific accretion onto two pre-existing dust populations
// (silicate + carbonaceous). Active when dust_species_track == 1.
//   - carbonaceous: rate-limited by gas-phase carbon
//   - silicate: rate-limited by min over {Mg, Si, Fe, O} of (rho_X / f_X),
//     following the Choban+2022 MNRAS 514, 4506 §2.2 key-reactant approach
// Per-species tau_accr uses the Hirashita 2011 section 2.6 normalization:
// n_H = 1e3 cm^-3, T = 50 K, S = 0.3, a = 0.1 micron, and solar key-species
// abundance, with S rescaled by dust_growth_sticking_coeff and a rescaled by
// dust_grainsize / 0.1. This is independent of the bulk SIMBA
// dust_growth_densref. No bulk dM, no partitioning — Phase D wires the species
// outputs into dust_update().
void dust_growth_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* growth_dM_silicate,  // output: silicate accretion rate
    double* growth_dM_carbon     // output: carbonaceous accretion rate
);

// Calculates dust destruction rates from SNe shocks and thermal sputtering.
// Stores the mass change dM for each cell in destruction_dM array.
void dust_destruction(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* t_gas,
    double* destruction_dM  // output: mass change rate for each cell
);

// Species-specific destruction (SN shocks + thermal sputtering) onto the two
// dust populations. Active when dust_species_track == 1.
//   - shock yield: graphite is baseline (factor 1.0); silicate follows the
//     Slavin+2015 standard SNR gas-cleared mass ratio 990/600 = 1.65
//     [REF: Slavin, Dwek, Jones 2015 ApJ 803, 7; Jones+1996 ApJ 469, 740]
//   - thermal sputtering: species-specific tau_ref
//     [REF (silicate): Tsai & Mathews 1995 ApJ 448, 84;
//      REF (carbon, ~2x silicate): Nozawa+2006 ApJ 648, 435]
//     with Draine & Salpeter 1979 / Tielens+1994 scaling form
// No bulk dM, no partitioning — Phase D wires the species outputs into
// dust_update().
void dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* t_gas,
    double* destruction_dM_silicate,  // output: silicate destruction rate
    double* destruction_dM_carbon     // output: carbonaceous destruction rate
);

// Update the density fields using calculated mass changes.
void dust_update(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value,
    const double* growth_dM,         // input: mass change from growth
    const double* destruction_dM,    // input: mass change from destruction
    bool dryrun);

// Species-specific field update for the two-species path (dust_species_track==1).
// Per-channel mass exchange:
//   - carbon channel:    rho_dust_carbonaceous <-> metal_density_carbon
//   - silicate channel:  rho_dust_silicate <-> {Mg, Fe, Si, O} at stoichiometric
//                        mass fractions f_X (Choban+2022 §2.2; 50/50 olivine +
//                        pyroxene mix, Draine 2003 / Dwek 1998).
// Per-channel pre-cap in absolute mass units replaces the legacy 3-way active[]
// shortfall mask. No SN injection here — Phase D drops the in-Grackle
// dust_creation pathway from the species branch; host code seeds dust species
// directly (or via inject_pathway machinery in make_consistent).
void dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value,
    const double* growth_dM_silicate,      // input: silicate accretion rate
    const double* growth_dM_carbon,        // input: carbonaceous accretion rate
    const double* destruction_dM_silicate, // input: silicate destruction rate (<=0)
    const double* destruction_dM_carbon,   // input: carbonaceous destruction rate (<=0)
    bool dryrun);

}  // namespace grackle::impl

#endif  // DUST_GROWTH_AND_DESTRUCTION_HPP
