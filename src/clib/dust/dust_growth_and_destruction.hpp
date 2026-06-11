#ifndef DUST_GROWTH_AND_DESTRUCTION_HPP
#define DUST_GROWTH_AND_DESTRUCTION_HPP

#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "support/index_helper.hpp"
#include "internal_units.hpp"
#include "phys_constants.h"
#include "fortran_func_decls.h"

namespace grackle::impl {

// Calculates dust growth rates (accretion) onto grain surfaces.
// Stores the mass change rate for each cell in growth_dM.
void dust_growth(chemistry_data* my_chemistry, grackle_field_data* my_fields,
                 InternalGrUnits internalu, IndexRange idx_range,
                 const gr_mask_type* itmask, const double* dt_value,
                 const double* t_gas,
                 double* growth_dM  // output: mass change rate for each cell
);

// Species-specific accretion onto three pre-existing dust populations
// (Mg-silicate + Fe-silicate + carbonaceous). Active when
// dust_species_track == 1.
//   - carbonaceous: bottlenecked by gas-phase carbon
//   - total silicate: bottlenecked by O, Si, and Mg+Fe, then split between
//     Mg2SiO4 and Fe2SiO4 endmembers according to local Mg/Fe abundance
// Uses the COLIBRE dense-gas accretion normalization with clumped
// n_H' = C n_H, T_ref = 10 K, n_H_ref = 10 cm^-3, S_ref = 0.3, a_ref = 0.1
// micron, and tau_ref values stored in Myr
// [REF: Trayford+2026 MNRAS 545, staf2040].
// The returned rates are equivalent rates whose product with dt gives the
// exponential e-folding mass update.
void dust_growth_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* growth_dM_mg_silicate,  // output: Mg-silicate accretion rate
    double* growth_dM_fe_silicate,  // output: Fe-silicate accretion rate
    double* growth_dM_carbon        // output: carbonaceous accretion rate
);

// Calculates dust destruction rates from SNe shocks and thermal sputtering.
// dt_full is the full external timestep (sne_rate is normalized per external
// step, not per subcycle). Stores the mass change rate for each cell in
// destruction_dM.
void dust_destruction(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, double dt_full, const double* t_gas,
    double* destruction_dM  // output: mass change rate for each cell
);

// Species-specific destruction (SN shocks + thermal sputtering) of the three
// dust populations. Active when dust_species_track == 1.
//   - shock yield: graphite is baseline (factor 1.0); both silicate endmembers
//     follow the Slavin+2015 standard SNR gas-cleared mass ratio
//     990/600 = 1.65
//     [REF: Slavin, Dwek, Jones 2015 ApJ 803, 7; Jones+1996 ApJ 469, 740]
//   - thermal sputtering: common COLIBRE/Tsai-Mathews tau_ref for all species
//     with (a/0.1) (n_H/cm^-3)^-1 [1 + (T/2e6 K)^-2.5] scaling
//     [REF: Tsai & Mathews 1995 ApJ 448, 84]
// Destruction uses exponential decay; the returned rates are equivalent rates
// whose product with dt gives the e-folding mass loss.
void dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, double dt_full, const double* t_gas,
    double* destruction_dM_mg_silicate,  // output: Mg-silicate destruction rate
    double* destruction_dM_fe_silicate,  // output: Fe-silicate destruction rate
    double* destruction_dM_carbon        // output: carbonaceous destruction rate
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
// (dust_species_track == 1). Per-channel mass exchange:
//   - carbon channel:    rho_dust_carbonaceous <-> metal_density_carbon
//   - Mg-sil channel:    rho_dust_mg_silicate <-> {Mg, Si, O} as Mg2SiO4
//   - Fe-sil channel:    rho_dust_fe_silicate <-> {Fe, Si, O} as Fe2SiO4
// Growth is capped by the limiting gas-phase reactant, destruction by the
// available dust mass. No SN injection here — host code seeds dust species
// directly (or via the inject_pathway machinery in make_consistent).
void dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value,
    const double* growth_dM_mg_silicate,  // input: Mg-silicate accretion rate
    const double* growth_dM_fe_silicate,  // input: Fe-silicate accretion rate
    const double* growth_dM_carbon,       // input: carbonaceous accretion rate
    const double* destruction_dM_mg_silicate,  // input: Mg-silicate destruction rate (<=0)
    const double* destruction_dM_fe_silicate,  // input: Fe-silicate destruction rate (<=0)
    const double* destruction_dM_carbon,       // input: carbonaceous destruction rate (<=0)
    bool dryrun);

}  // namespace grackle::impl

#endif  // DUST_GROWTH_AND_DESTRUCTION_HPP
