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

  // --- Element-resolved tracking setup ---
  bool track_elements = (my_chemistry->dust_model1_track_elements > 0 &&
                         my_fields->metal_density_carbon != nullptr &&
                         my_fields->metal_density_oxygen != nullptr &&
                         my_fields->dust_density_carbon != nullptr &&
                         my_fields->dust_density_oxygen != nullptr);

  grackle::impl::View<gr_float***> metal_C(
      track_elements ? my_fields->metal_density_carbon : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal_O(
      track_elements ? my_fields->metal_density_oxygen : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_C(
      track_elements ? my_fields->dust_density_carbon : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_O(
      track_elements ? my_fields->dust_density_oxygen : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust = dust(i, idx_range.j, idx_range.k);
      // metal_density holds the total gas-phase metals (C + O + others)
      double rho_metal_total = metal(i, idx_range.j, idx_range.k);
      double dt = dt_value[i];

      // Load element-resolved subsets (zero when not tracking)
      double rho_metal_C = 0.0, rho_metal_O = 0.0;
      double rho_dust_C = 0.0, rho_dust_O = 0.0;
      if (track_elements) {
        rho_metal_C = metal_C(i, idx_range.j, idx_range.k);
        rho_metal_O = metal_O(i, idx_range.j, idx_range.k);
        rho_dust_C  = dust_C(i, idx_range.j, idx_range.k);
        rho_dust_O  = dust_O(i, idx_range.j, idx_range.k);
      }
      // Non-C, non-O gas metals, derived from the total
      double rho_metal_other = rho_metal_total - rho_metal_C - rho_metal_O;

      double f_C     = (rho_metal_total > 0) ? rho_metal_C     / rho_metal_total : 0.0;
      double f_O     = (rho_metal_total > 0) ? rho_metal_O     / rho_metal_total : 0.0;
      double f_other = (rho_metal_total > 0) ? rho_metal_other / rho_metal_total : 0.0;
      double f_dust_C     = (rho_dust > 0) ? rho_dust_C / rho_dust : 0.0;
      double f_dust_O     = (rho_dust > 0) ? rho_dust_O / rho_dust : 0.0;
      double f_dust_other = (rho_dust > 0) ? 1.0 - (rho_dust_C + rho_dust_O) / rho_dust : 0.0;

      // --- Growth + destruction: exchanges mass between dust <-> metal ---
      // dM_exchange > 0 means net growth (metal -> dust)
      // dM_exchange < 0 means net destruction (dust -> metal)
      double dM_exchange = (growth_dM[i] + destruction_dM[i]) * dt;
      dM_exchange = std::max(-1.0 * rho_dust, dM_exchange);
      dM_exchange = std::min(0.9 * rho_metal_total, dM_exchange);

      // --- Stellar feedback injection (Dwek 1998 framework) ---
      double dM_create_dust = creation_dust_dM[i] * dt;
      double dM_create_metal = creation_metal_dM[i] * dt;

      // Apply conservation logic
      double dM_conserv = 0.0;
      if (rho_dust >= 0.0) {
        if (track_elements) {
        // --- GROWTH partitioning (Aoyama et al. 2017 proportional approach):
        // deplete gas-phase C, O, other metals proportionally to their

        // Compute intended subtractions from original dM_exchange
          double dM_C     = ((dM_exchange > 0) ? f_C : f_dust_C) * dM_exchange;
          double dM_O     = ((dM_exchange > 0) ? f_O : f_dust_O)     * dM_exchange;
          double dM_other = ((dM_exchange > 0) ? f_other : f_dust_other) * dM_exchange;

          // Cap each independently
          if (dM_C > ((dM_exchange > 0) ? rho_metal_C : -rho_dust_C)) {
              dM_C = (dM_exchange > 0) ? rho_metal_C : dM_C;
              rho_metal_C -= (dM_exchange > 0) ? rho_metal_C : dM_C;
          } else {
              rho_metal_C -= (dM_exchange > 0) ? dM_C : -rho_dust_C;
              dM_C = (dM_exchange > 0) ? dM_C : -rho_dust_C;
          }

          if (dM_O > ((dM_exchange > 0) ? rho_metal_O : -rho_dust_O)) {
              dM_O = (dM_exchange > 0) ? rho_metal_O : dM_O;
              rho_metal_O -= (dM_exchange > 0) ? rho_metal_O : dM_O;
          } else {
              rho_metal_O -= (dM_exchange > 0) ? dM_O : -rho_dust_O;
              dM_O = (dM_exchange > 0) ? dM_O : -rho_dust_O;
          }
          
          if (dM_other > ((dM_exchange > 0) ? rho_metal_other : -(rho_dust - rho_dust_C - rho_dust_O))) {
              dM_other = (dM_exchange > 0) ? rho_metal_other : dM_other;
              rho_metal_other -= (dM_exchange > 0) ? rho_metal_other : dM_other;
          } else {
              rho_metal_other -= (dM_exchange > 0) ? dM_other : -(rho_dust - rho_dust_C - rho_dust_O);
              dM_other = (dM_exchange > 0) ? dM_other : -(rho_dust - rho_dust_C - rho_dust_O);
          }

          // Dust receives exactly what was actually taken
          rho_dust_C     += dM_C;
          rho_dust_O     += dM_O;
          dM_exchange = dM_C + dM_O + dM_other;
          rho_dust = rho_dust + dM_exchange + dM_create_dust;
          rho_metal_C += dM_create_metal * f_C;
          rho_metal_O += dM_create_metal * f_O;
          rho_metal_other += dM_create_metal * f_other;

        } else {
          // Backward-compatible: all exchange with bulk metal_density
          rho_metal_other = rho_metal_other - dM_exchange + dM_create_metal;
          rho_dust = rho_dust + dM_exchange + dM_create_dust;

        }
      } else {
        dM_conserv = -rho_dust;
        rho_dust = 0.0;

        if (track_elements) {
          bool converged = false;
          if (dM_conserv > rho_metal_total) {
            // Not enough dust species to cover
            converged = false;
            fprintf(stderr,
                    "[DustConservation] WARNING: dust species deficit — "
                    "need %.6e but only %.6e available "
                    "(C=%.6e, O=%.6e). Capping transfer.\n",
                    dM_conserv, rho_metal_total,
                    rho_dust_C, rho_dust_O);
            std::exit(22);
          } else {
            if (rho_metal_C >= f_C * dM_conserv && rho_metal_O >= f_O * dM_conserv && rho_metal_other >= f_other * dM_conserv) {
              rho_metal_C -= f_C*dM_conserv;
              rho_metal_O -= f_O*dM_conserv;
              rho_metal_other -= f_other*dM_conserv;
            } else {
              bool active[3] = {true, true, true};
              double f_back[3] = {f_C, f_O, f_other};
              double* rho_metal_arr[3] = { &rho_metal_C, &rho_metal_O, &rho_metal_other };
              double dM_back[3]    = {f_C*dM_conserv,f_O*dM_conserv,f_other*dM_conserv};
              while (!converged) {
                converged = true;
                double shortfall = 0.0;
                double f_sum = 0.0;

                for (int j = 0; j < 3; j++) {
                  if (*rho_metal_arr[j] < dM_back[j]) {
                      shortfall += dM_back[j] - *rho_metal_arr[j];
                      dM_back[j] = *rho_metal_arr[j];
                      active[j] = false;
                  }
                }

                if (shortfall > 0.0) {
                  for (int j = 0; j < 3; j++) {
                    f_sum += active[j]*f_back[j];
                  }
                  if (f_sum > 0.0) {
                    converged = false;
                      for (int j = 0; j < 3; j++) {
                        if (active[j]) {
                          dM_back[j] += (f_back[j]/f_sum)*shortfall;
                        }
                      }
                  }
                }
              }
              rho_metal_C -= dM_back[0];
              rho_metal_O -= dM_back[1];
              rho_metal_other -= dM_back[2];
            }
          }
        } else {
            rho_metal_other -= dM_conserv;
        }
      }

      // Safety checks
      if (rho_dust < 0) {
        fprintf(stderr,
                "ERROR: Negative dust density at cell %d: rho_dust=%e\n", i,
                rho_dust);
        std::exit(21);
      }

      // Clamp element fields to prevent negative values from round-off
      if (track_elements) {
        rho_metal_C = std::max(0.0, rho_metal_C);
        rho_metal_O = std::max(0.0, rho_metal_O);
        rho_dust_C  = std::max(0.0, rho_dust_C);
        rho_dust_O  = std::max(0.0, rho_dust_O);
        // Ensure dust sub-components don't exceed total dust
        rho_dust_C = std::min(rho_dust_C, rho_dust);
        rho_dust_O = std::min(rho_dust_O, rho_dust - rho_dust_C);
      }

      // Gas density update:
      // density includes all metal fields as a subset, but NOT dust_density.
      // Track how total gas-phase metals changed to update gas density.
      // Compute after clamping so metal_C/O remain subsets of the total.
      double old_metal_total = metal(i, idx_range.j, idx_range.k);
      double new_metal_total = rho_metal_other + rho_metal_C + rho_metal_O;
      rho_gas = rho_gas + (new_metal_total - old_metal_total);

      fprintf(stderr,
              "internal: dt=%e growth=%.10e destruct=%.10e "
              "cre_dust=%.10e cre_metal=%.10e "
              "gas=%.15e dust=%.15e metal=%.15e\n",
              dt, growth_dM[i], destruction_dM[i],
              creation_dust_dM[i], creation_metal_dM[i],
              rho_gas, rho_dust, new_metal_total);
      // fprintf(stderr, "checking: %e\n", rho_gas + rho_dust);

      // Update the fields. metal_density now stores the *total* gas-phase
      // metals; metal_density_carbon / oxygen are subsets of it.
      if (!dryrun) {
        dust(i, idx_range.j, idx_range.k) = (gr_float)rho_dust;
        metal(i, idx_range.j, idx_range.k) = (gr_float)new_metal_total;
        d(i, idx_range.j, idx_range.k) = (gr_float)rho_gas;
        if (track_elements) {
          metal_C(i, idx_range.j, idx_range.k) = (gr_float)rho_metal_C;
          metal_O(i, idx_range.j, idx_range.k) = (gr_float)rho_metal_O;
          dust_C(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust_C;
          dust_O(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust_O;
        }
      }
    }
  }
}