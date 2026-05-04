#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "dust_growth_and_destruction.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

namespace {
const double sec_per_year = 3.155e7;

const double tiny_value = 1.0e-20;
const double huge_value = 1.0e+20;

const double t_ref = 20;

// Gate thresholds for dust evolution: skip cells where both dust and metals
// are negligible fractions of the baryon density.  The metal threshold matches
// the 1e-9 ratio used in make_consistent (make_consistent.cpp ~line 530) to
// distinguish metal-poor from metal-rich cells.  The dust threshold uses the
// same ratio for consistency — below 1 ppb of baryon density, any dust
// evolution would just integrate numerical noise.
const double dust_gate_threshold = 1.0e-9;
const double metal_gate_threshold = 1.0e-9;

}  // namespace

// ==========================================
// DUST GROWTH (ACCRETION)
// ==========================================
void grackle::impl::dust_growth(chemistry_data* my_chemistry,
                                grackle_field_data* my_fields,
                                InternalGrUnits internalu, IndexRange idx_range,
                                const gr_mask_type* itmask,
                                const double* dt_value, const double* t_gas,
                                double* growth_dM) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);
  double tau_ref =
      my_chemistry->dust_growth_tauref * 1e9 * sec_per_year / internalu.tbase1;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    // Initialize to zero
    growth_dM[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_dust = dust(i, idx_range.j, idx_range.k);
      // metal_density holds total gas-phase metals (C, O, and all others)
      double rho_metal = metal(i, idx_range.j, idx_range.k);
      double rho_d = d(i, idx_range.j, idx_range.k);

      // No metals to accrete onto grains — skip growth
      if (rho_metal < metal_gate_threshold * rho_d) {
        continue;
      }

      double temp = t_gas[i];
      double dt = dt_value[i];

      double tau_accr0 = tau_ref *
                         (my_chemistry->dust_growth_densref / dens_proper) *
                         std::pow(t_ref / temp, 0.5);
      double rho_metal_eff = std::max(rho_metal, tiny_value);
      double tau_accr =
          tau_accr0 * (my_chemistry->SolarMetalFractionByMass / rho_metal_eff);
      tau_accr = std::min(tau_accr, huge_value);
      tau_accr = std::max(tau_accr, tiny_value);
      double frac_metal_available = 0.0;
      if (rho_metal <= 0.0) {
        frac_metal_available = 0.0;
      } else if (rho_dust > 0.0 && rho_metal < 1e-12 * rho_dust) {
        frac_metal_available = rho_metal / rho_dust;
      } else {
        frac_metal_available = rho_metal / (rho_dust + rho_metal);
      }
      frac_metal_available = std::clamp(frac_metal_available, 0.0, 1.0);
      double growth_rate = frac_metal_available * (rho_dust / tau_accr);
      double dM = std::min(growth_rate, rho_metal / dt);

      // Store the calculated mass change in the output array
      growth_dM[i] = dM;
    }
  }
}

// ==========================================
// DUST GROWTH (SPECIES-SPECIFIC: silicate + carbonaceous)
// ==========================================
// Two parallel accretion rates onto pre-existing dust seeds, gated by
// dust_species_track == 1. Carbonaceous growth is rate-limited by gas-phase
// carbon. Silicate growth is rate-limited by the least-available key
// reactant in {Mg, Si, Fe, O} weighted by stoichiometric mass fraction
// f_X (Choban+2022 MNRAS 514, 4506 §2.2). The structural form of tau_accr
// follows Hirashita 2011 ApJ 743, 159 Eq. (16)-(17):
//   tau_accr = tau_ref · (rho_ref / rho) · (T_ref/T)^0.5 · (Z_sun / Z_X)
// reusing the same dust_growth_densref / SolarMetalFractionByMass scaling
// convention as the bulk dust_growth() above.
void grackle::impl::dust_growth_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* growth_dM_silicate, double* growth_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_sil(
      my_fields->dust_density_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_carb(
      my_fields->dust_density_carbonaceous, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mC(
      my_fields->metal_density_carbon, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mO(
      my_fields->metal_density_oxygen, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mMg(
      my_fields->metal_density_magnesium, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mSi(
      my_fields->metal_density_silicon, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mFe(
      my_fields->metal_density_iron, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);
  // tau_ref values are stored in Gyr (Hirashita 2011 Table 1 convention)
  double tau_ref_sil = my_chemistry->dust_growth_tauref_silicate * 1e9 *
                       sec_per_year / internalu.tbase1;
  double tau_ref_carb = my_chemistry->dust_growth_tauref_carbon * 1e9 *
                        sec_per_year / internalu.tbase1;

  double f_Mg = my_chemistry->dust_silicate_f_Mg;
  double f_Fe = my_chemistry->dust_silicate_f_Fe;
  double f_Si = my_chemistry->dust_silicate_f_Si;
  double f_O  = my_chemistry->dust_silicate_f_O;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    growth_dM_silicate[i] = 0.0;
    growth_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_sil_i  = dust_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);

      double rho_C  = mC(i, idx_range.j, idx_range.k);
      double rho_O  = mO(i, idx_range.j, idx_range.k);
      double rho_Mg = mMg(i, idx_range.j, idx_range.k);
      double rho_Si = mSi(i, idx_range.j, idx_range.k);
      double rho_Fe = mFe(i, idx_range.j, idx_range.k);

      double T  = t_gas[i];
      double dt = dt_value[i];

      double accr_struct = (my_chemistry->dust_growth_densref / dens_proper) *
                           std::pow(t_ref / T, 0.5);

      // ---------- Carbonaceous: rate-limited by gas-phase C ----------
      if (rho_C >= metal_gate_threshold * rho_gas && rho_dust_carb_i > 0.0) {
        double rho_C_eff = std::max(rho_C, tiny_value);
        double tau_accr_carb = tau_ref_carb * accr_struct *
                               (my_chemistry->SolarMetalFractionByMass /
                                rho_C_eff);
        tau_accr_carb = std::min(std::max(tau_accr_carb, tiny_value),
                                 huge_value);
        double frac_avail = rho_C / (rho_dust_carb_i + rho_C);
        frac_avail = std::clamp(frac_avail, 0.0, 1.0);
        double rate = frac_avail * (rho_dust_carb_i / tau_accr_carb);
        growth_dM_carbon[i] = std::min(rate, rho_C / dt);
      }

      // ---------- Silicate: Choban+2022 key-reactant ----------
      // For each key element X, the maximum silicate dust mass that could
      // be assembled from X is rho_X / f_X. The bottleneck element sets
      // the rate.
      double mass_from_Mg = (f_Mg > 0.0) ? rho_Mg / f_Mg : huge_value;
      double mass_from_Fe = (f_Fe > 0.0) ? rho_Fe / f_Fe : huge_value;
      double mass_from_Si = (f_Si > 0.0) ? rho_Si / f_Si : huge_value;
      double mass_from_O  = (f_O  > 0.0) ? rho_O  / f_O  : huge_value;
      double rho_sil_pool = std::min(std::min(mass_from_Mg, mass_from_Fe),
                                     std::min(mass_from_Si, mass_from_O));

      if (rho_sil_pool >= metal_gate_threshold * rho_gas &&
          rho_dust_sil_i > 0.0) {
        double rho_pool_eff = std::max(rho_sil_pool, tiny_value);
        double tau_accr_sil = tau_ref_sil * accr_struct *
                              (my_chemistry->SolarMetalFractionByMass /
                               rho_pool_eff);
        tau_accr_sil = std::min(std::max(tau_accr_sil, tiny_value),
                                huge_value);
        double frac_avail = rho_sil_pool / (rho_dust_sil_i + rho_sil_pool);
        frac_avail = std::clamp(frac_avail, 0.0, 1.0);
        double rate = frac_avail * (rho_dust_sil_i / tau_accr_sil);
        growth_dM_silicate[i] = std::min(rate, rho_sil_pool / dt);
      }
    }
  }
}

// ==========================================
// 2. DUST DESTRUCTION (SNe + SPUTTERING)
// ==========================================
void grackle::impl::dust_destruction(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* t_gas, double* destruction_dM) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  bool use_sne = (my_chemistry->use_sne_field > 0);
  grackle::impl::View<gr_float***> sne(
      use_sne ? my_fields->sne_rate : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  bool use_tau_dest = (my_chemistry->use_tau_dest_field > 0);
  grackle::impl::View<gr_float***> tau_dest_field(
      use_tau_dest ? my_fields->tau_dest : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);

  double Ms100 = 6800.0 * my_chemistry->sne_coeff *
                 (100.0 / my_chemistry->sne_shockspeed) *
                 (100.0 / my_chemistry->sne_shockspeed) * SolarMass /
                 (internalu.urho * std::pow(internalu.uxyz, 3));

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    // Initialize to zero
    destruction_dM[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust = dust(i, idx_range.j, idx_range.k);

      // No dust to destroy — skip destruction
      if (rho_dust < dust_gate_threshold * rho_gas) {
        continue;
      }

      double sne_this = use_sne ? sne(i, idx_range.j, idx_range.k) : 0.0;
      double temp = t_gas[i];
      double dt = dt_value[i];
      double tau_dest = 0;

      double dM = 0;
      double dM_shock = 0.0;

      if (use_tau_dest) {
        // user-provided destruction timescale
        tau_dest = tau_dest_field(i, idx_range.j, idx_range.k);
        if (tau_dest <= 0) {
          tau_dest = 1e20;
        }
        dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
      } else if (use_sne) {
        // destruction by SN shocks
        if (sne_this <= 0) {
          tau_dest = 1e20;
          // dM_shock = 0.0;
        } else {
          tau_dest = rho_gas /
                     (Ms100 * sne_this * my_chemistry->dust_destruction_eff) *
                     dt;
          dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
        }
      }

      if (temp >= std::pow(10, 5)) {
        // destruction by thermal sputtering
        double tau_sput = 1.7e8 * sec_per_year / internalu.tbase1 *
                          (my_chemistry->dust_grainsize / 0.1) *
                          (1.0e-27 / (dens_proper * rho_gas)) *
                          (std::pow((2.0e6 / temp), 2.5) + 1.0);

        if (dM_shock >= rho_dust / dt) {
          if (dM_shock > rho_dust / dt) {
            std::cout << "WARNING: dM_shock > M_dust SNe shock destruction, "
                      << sne_this << ", " << tau_dest << std::endl;
          }
        } else {
          dM_shock = dM_shock + rho_dust / tau_sput * 3.0;
          dM_shock = std::min(dM_shock, rho_dust / dt);
        }
      }
      // dM = - rho_dust * dM_shock;
      dM = -dM_shock;
      if (std::isnan(dM)) {
        std::cout << "dM calculated as NaN, " << dM << std::endl;
      }

      // Store the calculated mass change in the output array
      destruction_dM[i] = dM;
    }
  }
}

// ==========================================
// DUST DESTRUCTION (SPECIES-SPECIFIC: silicate + carbonaceous)
// ==========================================
// SN-shock destruction + thermal sputtering applied independently to each
// dust species, gated by dust_species_track == 1. Carbonaceous (graphite)
// is the shock-vulnerability baseline; silicates are ~1.5x more easily
// shattered by SN shocks under canonical ISM conditions
// [REF: Slavin, Dwek, Jones 2015 ApJ 803, 7; Jones+1996 ApJ 469, 740].
// Thermal sputtering uses species-specific tau_ref values
// [REF (silicate): Tsai & Mathews 1995 ApJ 448, 84;
//  REF (carbon, ~2x silicate): Nozawa+2006 ApJ 648, 435]
// with the same Draine & Salpeter 1979 / Tielens+1994 scaling form
// (a/0.1) * (1e-27/(rho*n)) * ((2e6/T)^2.5 + 1) used by the bulk path.
// Phase D wires the species outputs into dust_update().
void grackle::impl::dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* destruction_dM_silicate, double* destruction_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_sil(
      my_fields->dust_density_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_carb(
      my_fields->dust_density_carbonaceous, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  bool use_sne = (my_chemistry->use_sne_field > 0);
  grackle::impl::View<gr_float***> sne(
      use_sne ? my_fields->sne_rate : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  bool use_tau_dest = (my_chemistry->use_tau_dest_field > 0);
  grackle::impl::View<gr_float***> tau_dest_field(
      use_tau_dest ? my_fields->tau_dest : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);

  // Base shock-yield mass coefficient (Li+2019 framework).
  double Ms100 = 6800.0 * my_chemistry->sne_coeff *
                 (100.0 / my_chemistry->sne_shockspeed) *
                 (100.0 / my_chemistry->sne_shockspeed) * SolarMass /
                 (internalu.urho * std::pow(internalu.uxyz, 3));

  // Species-specific shock-vulnerability multipliers. Graphite is the
  // baseline (1.0); silicate is ~1.5x more easily destroyed under
  // canonical SN-shock conditions. [REF: Slavin+2015 ApJ 803, 7]
  const double shock_factor_carbon   = 1.0;
  const double shock_factor_silicate = 1.5;

  // Species-specific thermal sputtering tau_ref (params stored in years)
  double tau_sput_ref_sil = my_chemistry->dust_sputter_tauref_silicate *
                            sec_per_year / internalu.tbase1;
  double tau_sput_ref_carb = my_chemistry->dust_sputter_tauref_carbon *
                             sec_per_year / internalu.tbase1;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    destruction_dM_silicate[i] = 0.0;
    destruction_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_sil_i  = dust_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);

      bool sil_active  = (rho_dust_sil_i  >= dust_gate_threshold * rho_gas);
      bool carb_active = (rho_dust_carb_i >= dust_gate_threshold * rho_gas);
      if (!sil_active && !carb_active) continue;

      double sne_this = use_sne ? sne(i, idx_range.j, idx_range.k) : 0.0;
      double temp = t_gas[i];
      double dt = dt_value[i];

      // Common (species-independent) sputtering structural factor
      double sput_struct = (my_chemistry->dust_grainsize / 0.1) *
                           (1.0e-27 / (dens_proper * rho_gas)) *
                           (std::pow((2.0e6 / temp), 2.5) + 1.0);

      // Helper computing destruction rate for one species
      auto compute_dM = [&](double rho_dust, double shock_factor,
                            double tau_sput_ref) -> double {
        double dM_shock = 0.0;
        double tau_dest = 1e20;

        if (use_tau_dest) {
          tau_dest = tau_dest_field(i, idx_range.j, idx_range.k);
          if (tau_dest <= 0) tau_dest = 1e20;
          // Apply species multiplier on top of user-supplied tau_dest:
          // higher shock_factor -> shorter tau -> faster destruction.
          tau_dest = tau_dest / shock_factor;
          dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
        } else if (use_sne) {
          if (sne_this > 0.0) {
            tau_dest = rho_gas /
                       (Ms100 * shock_factor * sne_this *
                        my_chemistry->dust_destruction_eff) * dt;
            dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
          }
        }

        if (temp >= std::pow(10, 5)) {
          double tau_sput = tau_sput_ref * sput_struct;
          if (dM_shock < rho_dust / dt) {
            dM_shock = dM_shock + rho_dust / tau_sput * 3.0;
            dM_shock = std::min(dM_shock, rho_dust / dt);
          }
        }

        double dM = -dM_shock;
        if (std::isnan(dM)) {
          std::cout << "dM (species) calculated as NaN, " << dM << std::endl;
          dM = 0.0;
        }
        return dM;
      };

      if (carb_active) {
        destruction_dM_carbon[i] = compute_dM(
            rho_dust_carb_i, shock_factor_carbon, tau_sput_ref_carb);
      }
      if (sil_active) {
        destruction_dM_silicate[i] = compute_dM(
            rho_dust_sil_i, shock_factor_silicate, tau_sput_ref_sil);
      }
    }
  }
}

// ==========================================
// DUST UPDATE (SPECIES-SPECIFIC: silicate + carbonaceous)
// ==========================================
// Per-channel mass exchange between dust species and their corresponding
// gas-phase reactant pools.  Active when dust_species_track == 1.
//
//   carbon channel:   rho_dust_carbonaceous <-> rho_metal_carbon
//   silicate channel: rho_dust_silicate <-> {Mg, Fe, Si, O} at fractions f_X
//                     (Choban+2022 §2.2; 50/50 olivine+pyroxene per
//                      Draine 2003 / Dwek 1998).
//
// Per-channel pre-cap in absolute mass units replaces the legacy 3-way active[]
// shortfall iteration: growth is capped by the limiting reactant, destruction
// is capped by available dust mass.  No SN injection on this path — host code
// seeds the species directly.
void grackle::impl::dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* growth_dM_silicate,
    const double* growth_dM_carbon, const double* destruction_dM_silicate,
    const double* destruction_dM_carbon, bool dryrun) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_sil(
      my_fields->dust_density_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_carb(
      my_fields->dust_density_carbonaceous, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mC(
      my_fields->metal_density_carbon, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mO(
      my_fields->metal_density_oxygen, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mMg(
      my_fields->metal_density_magnesium, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mSi(
      my_fields->metal_density_silicon, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> mFe(
      my_fields->metal_density_iron, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  const double f_Mg = my_chemistry->dust_silicate_f_Mg;
  const double f_Fe = my_chemistry->dust_silicate_f_Fe;
  const double f_Si = my_chemistry->dust_silicate_f_Si;
  const double f_O  = my_chemistry->dust_silicate_f_O;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] == MASK_FALSE) continue;

    double rho_gas = d(i, idx_range.j, idx_range.k);
    double rho_dust_sil_i  = dust_sil(i, idx_range.j, idx_range.k);
    double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);
    double rho_metal_total = metal(i, idx_range.j, idx_range.k);
    double rho_C  = mC(i, idx_range.j, idx_range.k);
    double rho_O  = mO(i, idx_range.j, idx_range.k);
    double rho_Mg = mMg(i, idx_range.j, idx_range.k);
    double rho_Si = mSi(i, idx_range.j, idx_range.k);
    double rho_Fe = mFe(i, idx_range.j, idx_range.k);
    double dt = dt_value[i];

    // dM > 0 means net growth (gas reactants -> dust);
    // dM < 0 means net destruction (dust -> gas reactants).
    double dM_carb = (growth_dM_carbon[i] + destruction_dM_carbon[i]) * dt;
    double dM_sil  = (growth_dM_silicate[i] + destruction_dM_silicate[i]) * dt;

    // ---------- Per-channel pre-cap ----------
    // Carbon channel: bounded by rho_C (growth) or rho_dust_carb (destruction).
    if (dM_carb > 0.0) {
      dM_carb = std::min(dM_carb, rho_C);
    } else {
      dM_carb = std::max(dM_carb, -rho_dust_carb_i);
    }

    // Silicate channel: growth limited by least-available reactant per
    // stoichiometric coefficient; destruction limited by silicate dust mass.
    if (dM_sil > 0.0) {
      double cap = dM_sil;
      if (f_Mg > 0.0) cap = std::min(cap, rho_Mg / f_Mg);
      if (f_Fe > 0.0) cap = std::min(cap, rho_Fe / f_Fe);
      if (f_Si > 0.0) cap = std::min(cap, rho_Si / f_Si);
      if (f_O  > 0.0) cap = std::min(cap, rho_O  / f_O);
      dM_sil = std::max(0.0, cap);
    } else {
      dM_sil = std::max(dM_sil, -rho_dust_sil_i);
    }

    // ---------- Apply ----------
    rho_dust_carb_i += dM_carb;
    rho_C           -= dM_carb;

    rho_dust_sil_i += dM_sil;
    rho_Mg -= dM_sil * f_Mg;
    rho_Fe -= dM_sil * f_Fe;
    rho_Si -= dM_sil * f_Si;
    rho_O  -= dM_sil * f_O;

    // Bulk dust = silicate + carbonaceous (Phase E will enforce as invariant
    // in make_consistent; computing it here keeps the bulk field in sync).
    double rho_dust_new = rho_dust_sil_i + rho_dust_carb_i;

    // Total gas-phase metal change: -(dM_carb + dM_sil) (metals lost when
    // dust grows).  metal_density_other (= total - C - O - Mg - Si - Fe) is
    // unchanged on the dust loop per Choban+2022 §2.2.
    double delta_metal_total = -(dM_carb + dM_sil);
    double rho_metal_new = rho_metal_total + delta_metal_total;

    // Gas density tracks metals as a subset (dust is not part of `density`).
    rho_gas += delta_metal_total;

    // Floors / NaN guard
    rho_dust_carb_i = std::max(0.0, rho_dust_carb_i);
    rho_dust_sil_i  = std::max(0.0, rho_dust_sil_i);
    rho_dust_new    = std::max(0.0, rho_dust_new);
    rho_C  = std::max(0.0, rho_C);
    rho_O  = std::max(0.0, rho_O);
    rho_Mg = std::max(0.0, rho_Mg);
    rho_Si = std::max(0.0, rho_Si);
    rho_Fe = std::max(0.0, rho_Fe);
    rho_metal_new = std::max(0.0, rho_metal_new);

    if (std::isnan(rho_dust_new) || std::isnan(rho_metal_new)) {
      std::cout << "dust_update_species: NaN at cell " << i
                << " dM_carb=" << dM_carb << " dM_sil=" << dM_sil << std::endl;
      continue;
    }

    if (!dryrun) {
      dust(i, idx_range.j, idx_range.k)      = (gr_float)rho_dust_new;
      dust_sil(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust_sil_i;
      dust_carb(i, idx_range.j, idx_range.k) = (gr_float)rho_dust_carb_i;
      metal(i, idx_range.j, idx_range.k)     = (gr_float)rho_metal_new;
      mC(i, idx_range.j, idx_range.k)        = (gr_float)rho_C;
      mO(i, idx_range.j, idx_range.k)        = (gr_float)rho_O;
      mMg(i, idx_range.j, idx_range.k)       = (gr_float)rho_Mg;
      mSi(i, idx_range.j, idx_range.k)       = (gr_float)rho_Si;
      mFe(i, idx_range.j, idx_range.k)       = (gr_float)rho_Fe;
      d(i, idx_range.j, idx_range.k)         = (gr_float)rho_gas;
    }
  }
}

void grackle::impl::dust_update(chemistry_data* my_chemistry,
                                grackle_field_data* my_fields,
                                InternalGrUnits internalu, IndexRange idx_range,
                                const gr_mask_type* itmask,
                                const double* dt_value, const double* growth_dM,
                                const double* destruction_dM, bool dryrun) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] == MASK_FALSE) continue;

    double rho_gas = d(i, idx_range.j, idx_range.k);
    double rho_dust = dust(i, idx_range.j, idx_range.k);
    double rho_metal_total = metal(i, idx_range.j, idx_range.k);
    double dt = dt_value[i];

    // dM_exchange > 0: net growth (metal -> dust); < 0: destruction (dust -> metal)
    double dM_exchange = (growth_dM[i] + destruction_dM[i]) * dt;
    dM_exchange = std::max(-1.0 * rho_dust, dM_exchange);
    dM_exchange = std::min(0.9 * rho_metal_total, dM_exchange);

    rho_dust         += dM_exchange;
    rho_metal_total  -= dM_exchange;
    rho_gas          -= dM_exchange;  // gas tracks metals as a subset

    rho_dust        = std::max(0.0, rho_dust);
    rho_metal_total = std::max(0.0, rho_metal_total);

    if (!dryrun) {
      dust(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust;
      metal(i, idx_range.j, idx_range.k) = (gr_float)rho_metal_total;
      d(i, idx_range.j, idx_range.k)     = (gr_float)rho_gas;
    }
  }
}