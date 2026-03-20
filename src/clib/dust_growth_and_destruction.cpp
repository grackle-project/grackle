#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "dust_growth_and_destruction.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

namespace {
const double k_boltz = 1.3806504e-16;
const double m_proton = 1.67262171e-24;
const double pi_val = 3.141592653589793;
const double sec_per_year = 3.155e7;

const double tiny_value = 1.0e-20;
const double huge_value = 1.0e+20;

const double t_ref = 20;

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
      double rho_metal = metal(i, idx_range.j, idx_range.k);
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
// DUST CREATION (STELLAR FEEDBACK)
// ==========================================
// Following the standard dust evolution framework (Dwek 1998, ApJ 501, 643;
// reviewed in Galliano et al. 2018, ARA&A 56, 673). Each SN injects a total
// metal mass m_Z (sne_metal_yield). A fraction delta (dust_condensation_eff)
// condenses into dust grains, while the remaining (1 - delta) fraction enters
// the gas phase as metals:
//
//   dM_dust/dt   = delta * m_Z * R_SN           (Eq. 1 of review 2504.10585)
//   dM_metal/dt  = (1 - delta) * m_Z * R_SN
//
// where R_SN = sne_rate (SN rate per unit volume per unit time).
//
void grackle::impl::dust_creation(chemistry_data* my_chemistry,
                                  grackle_field_data* my_fields,
                                  InternalGrUnits internalu, IndexRange idx_range,
                                  const gr_mask_type* itmask,
                                  const double* dt_value,
                                  double* creation_dust_dM,
                                  double* creation_metal_dM) {
  bool use_sne = (my_chemistry->use_sne_field > 0);
  if (!use_sne) {
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
      creation_dust_dM[i] = 0.0;
      creation_metal_dM[i] = 0.0;
    }
    return;
  }

  grackle::impl::View<gr_float***> sne(
      my_fields->sne_rate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Metal yield per SN in code mass units
  double M_yield_code = my_chemistry->sne_metal_yield * SolarMass /
                        (internalu.urho * std::pow(internalu.uxyz, 3));

  double f_cond = my_chemistry->dust_condensation_eff;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    creation_dust_dM[i] = 0.0;
    creation_metal_dM[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double sne_this = sne(i, idx_range.j, idx_range.k);
      double dt = dt_value[i];

      if (sne_this > 0.0) {
        double total_metal_inject = M_yield_code * sne_this / dt;
        creation_dust_dM[i] = f_cond * total_metal_inject;
        creation_metal_dM[i] = (1.0 - f_cond) * total_metal_inject;
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

void grackle::impl::dust_update(chemistry_data* my_chemistry,
                                grackle_field_data* my_fields,
                                InternalGrUnits internalu, IndexRange idx_range,
                                const gr_mask_type* itmask,
                                const double* dt_value, const double* growth_dM,
                                const double* destruction_dM,
                                const double* creation_dust_dM,
                                const double* creation_metal_dM,
                                bool dryrun) {
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
    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust = dust(i, idx_range.j, idx_range.k);
      double rho_metal = metal(i, idx_range.j, idx_range.k);
      double dt = dt_value[i];

      // --- Growth + destruction: exchanges mass between dust <-> metal ---
      double dM_exchange = (growth_dM[i] + destruction_dM[i]) * dt;
      dM_exchange = std::max(-1.0 * rho_dust, dM_exchange);
      dM_exchange = std::min(0.9 * rho_metal, dM_exchange);

      // --- Stellar feedback injection (Dwek 1998 framework) ---
      // SN ejecta inject total metal mass m_Z * R_SN * dt:
      //   delta * m_Z -> dust     (separate from gas density)
      //   (1-delta) * m_Z -> gas-phase metals (subset of gas density)
      double dM_create_dust = creation_dust_dM[i] * dt;
      double dM_create_metal = creation_metal_dM[i] * dt;

      // Apply conservation logic (from original code)
      double dM_conserv = 0.0;
      if (rho_dust >= 0.0) {
        rho_dust = rho_dust + dM_exchange + dM_create_dust;
        rho_metal = rho_metal - dM_exchange + dM_create_metal;
      } else {
        dM_conserv = rho_dust;
        rho_dust = rho_dust - dM_conserv;
        rho_metal = rho_metal + dM_conserv;
      }

      // Gas density update:
      // density includes metal_density as a subset, but NOT dust_density.
      // Exchange: metal<->dust moves mass in/out of gas, so density
      //   changes by the same amount as metal (-dM_exchange).
      // Stellar injection: only (1-delta)*m_Z enters gas (as metals).
      //   The delta*m_Z portion goes directly to dust (outside gas).
      // Both effects are captured by tracking how metal changed:
      double old_metal = metal(i, idx_range.j, idx_range.k);
      rho_gas = rho_gas + (rho_metal - old_metal);

      // Safety checks
      if (rho_dust < 0) {
        fprintf(stderr,
                "ERROR: Negative dust density at cell %d: rho_dust=%e\n", i,
                rho_dust);
        std::exit(21);
      }

      fprintf(stderr,
              "internal: dt=%e growth=%.10e destruct=%.10e "
              "cre_dust=%.10e cre_metal=%.10e "
              "gas=%.15e dust=%.15e metal=%.15e\n",
              dt, growth_dM[i], destruction_dM[i],
              creation_dust_dM[i], creation_metal_dM[i],
              rho_gas, rho_dust, rho_metal);

      // Update the fields
      if (!dryrun) {
        dust(i, idx_range.j, idx_range.k) = (gr_float)rho_dust;
        metal(i, idx_range.j, idx_range.k) = (gr_float)rho_metal;
        d(i, idx_range.j, idx_range.k) = (gr_float)rho_gas;
      }
    }
  }
}