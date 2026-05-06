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

// Species-specific accretion onto three pre-existing dust populations
// (olivine + pyroxene + carbonaceous). Active when dust_species_track == 1.
//   - carbonaceous: rate-limited by gas-phase carbon
//   - olivine: rate-limited by min over {Mg, Fe, Si, O} of (rho_X / f_X)
//   - pyroxene: rate-limited by min over {Mg, Si, O} of (rho_X / f_X)
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
    double* growth_dM_olivine,   // output: olivine accretion rate
    double* growth_dM_pyroxene,  // output: pyroxene accretion rate
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

// Species-specific destruction (SN shocks + thermal sputtering) onto the three
// dust populations. Active when dust_species_track == 1.
//   - shock yield: graphite is baseline (factor 1.0); olivine and pyroxene
//     both follow the
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
    double* destruction_dM_olivine,   // output: olivine destruction rate
    double* destruction_dM_pyroxene,  // output: pyroxene destruction rate
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

// Species-specific field update for the split-silicate path
// (dust_species_track==1).
// Per-channel mass exchange:
//   - carbon channel:    rho_dust_carbonaceous <-> metal_density_carbon
//   - olivine channel:   rho_dust_olivine <-> {Mg, Fe, Si, O} as MgFeSiO4
//   - pyroxene channel:  rho_dust_pyroxene <-> {Mg, Si, O} as MgSiO3
// Per-channel pre-cap in absolute mass units replaces the legacy 3-way active[]
// shortfall mask. No SN injection here — Phase D drops the in-Grackle
// dust_creation pathway from the species branch; host code seeds dust species
// directly (or via inject_pathway machinery in make_consistent).
void dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value,
    const double* growth_dM_olivine,       // input: olivine accretion rate
    const double* growth_dM_pyroxene,      // input: pyroxene accretion rate
    const double* growth_dM_carbon,        // input: carbonaceous accretion rate
    const double* destruction_dM_olivine,  // input: olivine destruction rate (<=0)
    const double* destruction_dM_pyroxene, // input: pyroxene destruction rate (<=0)
    const double* destruction_dM_carbon,   // input: carbonaceous destruction rate (<=0)
    bool dryrun);

}  // namespace grackle::impl

#endif  // DUST_GROWTH_AND_DESTRUCTION_HPP
