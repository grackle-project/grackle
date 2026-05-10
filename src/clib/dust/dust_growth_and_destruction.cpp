#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "dust_growth_and_destruction.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

namespace {
const double sec_per_year = 3.155e7;
const double sec_per_Myr_local = 1.0e6 * sec_per_year;

const double tiny_value = 1.0e-20;
const double huge_value = 1.0e+20;

const double t_ref = 20;
const double colibre_growth_t_ref = 10.0;
const double colibre_growth_nH_ref = 10.0;
const double colibre_sputter_t_ref = 2.0e6;

const double atomic_C  = 12.01;
const double atomic_O  = 16.00;
const double atomic_Mg = 24.31;
const double atomic_Si = 28.09;
const double atomic_Fe = 55.85;

double limited_expm1(double x) {
  if (!std::isfinite(x) || x <= 0.0) return 0.0;
  return std::expm1(std::min(x, 50.0));
}

double e_fold_growth_rate(double rho_dust, double dt, double tau) {
  if (rho_dust <= 0.0 || dt <= 0.0 || tau <= 0.0) return 0.0;
  return rho_dust * limited_expm1(dt / tau) / dt;
}

double e_fold_loss_rate(double rho_dust, double dt, double inv_tau) {
  if (rho_dust <= 0.0 || dt <= 0.0 || inv_tau <= 0.0) return 0.0;
  double x = std::min(dt * inv_tau, 50.0);
  return rho_dust * (-std::expm1(-x)) / dt;
}

// Solar metal mass fractions for the tracked dust-forming elements, relative
// to total solar metals. These match gracklepy.utilities.convenience, which
// seeds the dust_species_track gas reservoirs from metal_density.
const double solar_frac_C  = 0.15925782394660776;
const double solar_frac_O  = 0.4242932765702842;
const double solar_frac_Mg = 0.045644817372018066;
const double solar_frac_Si = 0.052744600629574714;
const double solar_frac_Fe = 0.08523143041944482;

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
// DUST GROWTH (SPECIES-SPECIFIC: Mg-silicate + Fe-silicate + carbonaceous)
// ==========================================
// COLIBRE-style accretion (Trayford+2026 section 3.3). Carbonaceous growth is
// bottlenecked by gas-phase C. Silicate growth is computed once for the total
// silicate reservoir using the maximum depletion factor over O, Si, and Mg+Fe,
// then split between Mg2SiO4 and Fe2SiO4 endmembers according to the local
// gas-phase Mg/Fe number abundance. The dense-gas rate uses a clumped
// hydrogen density n_H' = C n_H; thermal sputtering below deliberately uses the
// unclumped n_H. This species-tracking branch applies growth as an exponential
// e-folding update, returning an equivalent per-time rate for dust_update().
void grackle::impl::dust_growth_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* growth_dM_mg_silicate, double* growth_dM_fe_silicate,
    double* growth_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_mg_sil(
      my_fields->dust_density_mg_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_fe_sil(
      my_fields->dust_density_fe_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_carb(
      my_fields->dust_density_carbonaceous, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  bool use_H_fields = (my_chemistry->primordial_chemistry > 0);
  bool use_H2_fields = (my_chemistry->primordial_chemistry > 1);
  bool use_HD_fields = (my_chemistry->primordial_chemistry > 2);
  bool use_HeH_fields = (my_chemistry->primordial_chemistry > 3);
  grackle::impl::View<gr_float***> HI(
      use_H_fields ? my_fields->HI_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      use_H_fields ? my_fields->HII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(
      use_H2_fields ? my_fields->HM_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(
      use_H2_fields ? my_fields->H2I_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      use_H2_fields ? my_fields->H2II_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(
      use_HD_fields ? my_fields->HDI_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(
      use_HeH_fields ? my_fields->HDII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(
      use_HeH_fields ? my_fields->HeHII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
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
  // tau_ref values are stored in Myr (COLIBRE Table 2).
  double tau_ref_sil = my_chemistry->dust_growth_tauref_silicate *
                       sec_per_Myr_local / internalu.tbase1;
  double tau_ref_carb = my_chemistry->dust_growth_tauref_carbon *
                        sec_per_Myr_local / internalu.tbase1;
  double sticking_factor =
      (my_chemistry->dust_growth_sticking_coeff > 0.0)
          ? 0.3 / my_chemistry->dust_growth_sticking_coeff
          : huge_value;
  double grain_size_factor =
      (my_chemistry->dust_grainsize > 0.0)
          ? my_chemistry->dust_grainsize / 0.1
          : huge_value;

  double solar_nonmetal =
      std::max(1.0 - my_chemistry->SolarMetalFractionByMass, tiny_value);
  double solar_H = std::max(my_chemistry->HydrogenFractionByMass *
                            solar_nonmetal, tiny_value);
  double solar_mass_C  = my_chemistry->SolarMetalFractionByMass *
                         solar_frac_C / solar_H;
  double solar_mass_O  = my_chemistry->SolarMetalFractionByMass *
                         solar_frac_O / solar_H;
  double solar_mass_Mg = my_chemistry->SolarMetalFractionByMass *
                         solar_frac_Mg / solar_H;
  double solar_mass_Si = my_chemistry->SolarMetalFractionByMass *
                         solar_frac_Si / solar_H;
  double solar_mass_Fe = my_chemistry->SolarMetalFractionByMass *
                         solar_frac_Fe / solar_H;
  double solar_eps_C = solar_mass_C / atomic_C;
  double solar_eps_O = solar_mass_O / atomic_O;
  double solar_eps_Mg = solar_mass_Mg / atomic_Mg;
  double solar_eps_Si = solar_mass_Si / atomic_Si;
  double solar_eps_Fe = solar_mass_Fe / atomic_Fe;
  double solar_eps_MgFe = solar_eps_Mg + solar_eps_Fe;

  auto clumping_factor = [&](double nH) -> double {
    double cmax = std::max(my_chemistry->dust_growth_clumping_factor_max, 1.0);
    double nmin = std::max(my_chemistry->dust_growth_clumping_nH_min,
                           tiny_value);
    double nmax = std::max(my_chemistry->dust_growth_clumping_nH_max,
                           nmin * (1.0 + 1.0e-12));
    if (nH <= nmin || cmax <= 1.0) return 1.0;
    if (nH >= nmax) return cmax;
    double x = (std::log10(nH) - std::log10(nmin)) /
               (std::log10(nmax) - std::log10(nmin));
    return std::pow(cmax, std::clamp(x, 0.0, 1.0));
  };

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    growth_dM_mg_silicate[i] = 0.0;
    growth_dM_fe_silicate[i] = 0.0;
    growth_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_mg_sil_i = dust_mg_sil(i, idx_range.j, idx_range.k);
      double rho_dust_fe_sil_i = dust_fe_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);
      double rho_metal_total = metal(i, idx_range.j, idx_range.k);

      double rho_C  = mC(i, idx_range.j, idx_range.k);
      double rho_O  = mO(i, idx_range.j, idx_range.k);
      double rho_Mg = mMg(i, idx_range.j, idx_range.k);
      double rho_Si = mSi(i, idx_range.j, idx_range.k);
      double rho_Fe = mFe(i, idx_range.j, idx_range.k);

      double T  = std::max(t_gas[i], tiny_value);
      double dt = dt_value[i];
      double rho_nonmetal = rho_gas - rho_metal_total;
      if (rho_nonmetal <= 0.0) {
        continue;
      }
      double rho_H_nuclei = my_chemistry->HydrogenFractionByMass *
                            rho_nonmetal;
      if (use_H_fields) {
        rho_H_nuclei = HI(i, idx_range.j, idx_range.k) +
                       HII(i, idx_range.j, idx_range.k);
        if (use_H2_fields) {
          rho_H_nuclei += HM(i, idx_range.j, idx_range.k) +
                          H2I(i, idx_range.j, idx_range.k) +
                          H2II(i, idx_range.j, idx_range.k);
        }
        if (use_HD_fields) {
          rho_H_nuclei += HDI(i, idx_range.j, idx_range.k) / 3.0;
        }
        if (use_HeH_fields) {
          rho_H_nuclei += HDII(i, idx_range.j, idx_range.k) / 3.0 +
                          HeHII(i, idx_range.j, idx_range.k) / 5.0;
        }
      }
      double nH = rho_H_nuclei * dens_proper / mh;
      if (rho_gas <= 0.0 || rho_H_nuclei <= 0.0 || nH <= 0.0) {
        continue;
      }

      double nH_eff = nH * clumping_factor(nH);
      double accr_struct = (colibre_growth_nH_ref / nH_eff) *
                           std::pow(colibre_growth_t_ref / T, 0.5) *
                           grain_size_factor * sticking_factor;
      double inv_rho_H = 1.0 / rho_H_nuclei;
      double eps_C  = rho_C  * inv_rho_H / atomic_C;
      double eps_O  = rho_O  * inv_rho_H / atomic_O;
      double eps_Mg = rho_Mg * inv_rho_H / atomic_Mg;
      double eps_Si = rho_Si * inv_rho_H / atomic_Si;
      double eps_Fe = rho_Fe * inv_rho_H / atomic_Fe;
      double eps_MgFe = eps_Mg + eps_Fe;

      auto abundance_ratio = [](double solar_eps, double local_eps) {
        if (solar_eps <= 0.0) return 0.0;
        if (local_eps <= 0.0) return huge_value;
        return solar_eps / local_eps;
      };

      // ---------- Carbonaceous: rate-limited by gas-phase C ----------
      if (rho_C >= metal_gate_threshold * rho_gas &&
          rho_dust_carb_i > 0.0 && solar_eps_C > 0.0 && eps_C > 0.0 &&
          dt > 0.0) {
        double tau_accr_carb = tau_ref_carb * accr_struct *
                               abundance_ratio(solar_eps_C, eps_C);
        tau_accr_carb = std::min(std::max(tau_accr_carb, tiny_value),
                                 huge_value);
        double rate = e_fold_growth_rate(rho_dust_carb_i, dt,
                                         tau_accr_carb);
        growth_dM_carbon[i] = std::min(rate, rho_C / dt);
      }

      // ---------- Silicate: COLIBRE Mg+Fe composite bottleneck ----------
      double rho_dust_sil_i = rho_dust_mg_sil_i + rho_dust_fe_sil_i;
      if (rho_dust_sil_i > 0.0 && dt > 0.0 &&
          eps_O > 0.0 && eps_Si > 0.0 && eps_MgFe > 0.0) {
        double sil_ratio = std::max({
            abundance_ratio(solar_eps_O, eps_O),
            abundance_ratio(solar_eps_Si, eps_Si),
            abundance_ratio(solar_eps_MgFe, eps_MgFe)});
        double tau_accr = tau_ref_sil * accr_struct *
                          sil_ratio;
        tau_accr = std::min(std::max(tau_accr, tiny_value), huge_value);
        double sil_rate = e_fold_growth_rate(rho_dust_sil_i, dt,
                                             tau_accr);

        // Split total silicate accretion by endmember molecule abundance.
        // Mass fractions include the different Mg2SiO4 / Fe2SiO4 molecule
        // weights, so equal Mg and Fe number abundance reproduces the
        // COLIBRE equal-molecule seed split.
        double mg_weight = eps_Mg * (
            2.0 * atomic_Mg + atomic_Si + 4.0 * atomic_O);
        double fe_weight = eps_Fe * (
            2.0 * atomic_Fe + atomic_Si + 4.0 * atomic_O);
        double endmember_weight = mg_weight + fe_weight;
        if (endmember_weight > 0.0) {
          double mg_frac = mg_weight / endmember_weight;
          mg_frac = std::clamp(mg_frac, 0.0, 1.0);
          growth_dM_mg_silicate[i] = sil_rate * mg_frac;
          growth_dM_fe_silicate[i] = sil_rate * (1.0 - mg_frac);
        }
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
// DUST DESTRUCTION (SPECIES-SPECIFIC: Mg-silicate + Fe-silicate + carbonaceous)
// ==========================================
// SN-shock destruction + thermal sputtering applied independently to each
// dust species, gated by dust_species_track == 1. Carbonaceous (graphite) is
// the shock-vulnerability baseline; both silicate endmembers use the same
// silicate shock coefficient. Slavin+2015 Table 2 gives gas-cleared masses of
// 990 Msun for silicates and 600 Msun for carbonaceous grains in their
// standard SNR model, so silicates are destroyed about 1.65x faster
// [REF: Slavin, Dwek, Jones 2015 ApJ 803, 7; Jones+1996 ApJ 469, 740].
// Thermal sputtering follows the species-independent COLIBRE/Tsai-Mathews
// timescale tau_sp = 0.85 Myr (a/0.1) (n_H/cm^-3)^-1
// [1 + (T/2e6 K)^-2.5]. Species destruction is applied as exponential decay,
// returning an equivalent per-time rate for dust_update().
// Phase D wires the species outputs into dust_update().
void grackle::impl::dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* destruction_dM_mg_silicate, double* destruction_dM_fe_silicate,
    double* destruction_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_mg_sil(
      my_fields->dust_density_mg_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_fe_sil(
      my_fields->dust_density_fe_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_carb(
      my_fields->dust_density_carbonaceous, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  bool use_H_fields = (my_chemistry->primordial_chemistry > 0);
  bool use_H2_fields = (my_chemistry->primordial_chemistry > 1);
  bool use_HD_fields = (my_chemistry->primordial_chemistry > 2);
  bool use_HeH_fields = (my_chemistry->primordial_chemistry > 3);
  grackle::impl::View<gr_float***> HI(
      use_H_fields ? my_fields->HI_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HII(
      use_H_fields ? my_fields->HII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HM(
      use_H2_fields ? my_fields->HM_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2I(
      use_H2_fields ? my_fields->H2I_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> H2II(
      use_H2_fields ? my_fields->H2II_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDI(
      use_HD_fields ? my_fields->HDI_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HDII(
      use_HeH_fields ? my_fields->HDII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> HeHII(
      use_HeH_fields ? my_fields->HeHII_density : my_fields->density,
      my_fields->grid_dimension[0], my_fields->grid_dimension[1],
      my_fields->grid_dimension[2]);
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
  // baseline (1.0); silicate follows the standard Slavin+2015 SNR
  // gas-cleared mass ratio, 990/600 = 1.65.
  const double shock_factor_carbon   = 1.0;
  const double shock_factor_silicate = 1.65;

  // COLIBRE/Tsai-Mathews thermal sputtering tau_ref (stored in Myr).
  double tau_sput_ref = my_chemistry->dust_sputter_tauref *
                        sec_per_Myr_local / internalu.tbase1;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    destruction_dM_mg_silicate[i] = 0.0;
    destruction_dM_fe_silicate[i] = 0.0;
    destruction_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_mg_sil_i = dust_mg_sil(i, idx_range.j, idx_range.k);
      double rho_dust_fe_sil_i = dust_fe_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);
      double rho_metal_total = metal(i, idx_range.j, idx_range.k);

      bool mg_sil_active =
          (rho_dust_mg_sil_i >= dust_gate_threshold * rho_gas);
      bool fe_sil_active =
          (rho_dust_fe_sil_i >= dust_gate_threshold * rho_gas);
      bool carb_active = (rho_dust_carb_i >= dust_gate_threshold * rho_gas);
      if (!mg_sil_active && !fe_sil_active && !carb_active) continue;

      double sne_this = use_sne ? sne(i, idx_range.j, idx_range.k) : 0.0;
      if (rho_gas <= 0.0) continue;

      double temp = std::max(t_gas[i], tiny_value);
      double dt = dt_value[i];
      if (dt <= 0.0) continue;

      double rho_nonmetal = rho_gas - rho_metal_total;
      if (rho_nonmetal <= 0.0) continue;
      double rho_H_nuclei = my_chemistry->HydrogenFractionByMass *
                            rho_nonmetal;
      if (use_H_fields) {
        rho_H_nuclei = HI(i, idx_range.j, idx_range.k) +
                       HII(i, idx_range.j, idx_range.k);
        if (use_H2_fields) {
          rho_H_nuclei += HM(i, idx_range.j, idx_range.k) +
                          H2I(i, idx_range.j, idx_range.k) +
                          H2II(i, idx_range.j, idx_range.k);
        }
        if (use_HD_fields) {
          rho_H_nuclei += HDI(i, idx_range.j, idx_range.k) / 3.0;
        }
        if (use_HeH_fields) {
          rho_H_nuclei += HDII(i, idx_range.j, idx_range.k) / 3.0 +
                          HeHII(i, idx_range.j, idx_range.k) / 5.0;
        }
      }
      double nH = rho_H_nuclei * dens_proper / mh;
      if (rho_H_nuclei <= 0.0 || nH <= 0.0) continue;

      // Common (species-independent) sputtering structural factor.
      double sput_struct = (my_chemistry->dust_grainsize / 0.1) *
                           (1.0 / nH) *
                           (1.0 + std::pow(temp / colibre_sputter_t_ref,
                                           -2.5));

      // Helper computing the equivalent per-time destruction rate for one
      // species. Its dt product is an exponential e-folding mass loss.
      auto compute_dM = [&](double rho_dust, double shock_factor) -> double {
        double inv_tau_loss = 0.0;
        double tau_dest = 1e20;

        if (use_tau_dest) {
          tau_dest = tau_dest_field(i, idx_range.j, idx_range.k);
          if (tau_dest <= 0) tau_dest = 1e20;
          // Apply species multiplier on top of user-supplied tau_dest:
          // higher shock_factor -> shorter tau -> faster destruction.
          tau_dest = tau_dest / shock_factor;
          if (tau_dest > 0.0 && std::isfinite(tau_dest)) {
            inv_tau_loss += 1.0 / tau_dest;
          }
        } else if (use_sne) {
          if (sne_this > 0.0) {
            // sne_rate is a rate per code time, so add the physical inverse
            // shock-destruction timescale here. The exponential update below
            // supplies the dt dependence.
            double inv_tau_shock =
                Ms100 * shock_factor * sne_this *
                my_chemistry->dust_destruction_eff / (rho_gas*dt);
            if (inv_tau_shock > 0.0 && std::isfinite(inv_tau_shock)) {
              inv_tau_loss += inv_tau_shock;
            }
          }
        }

        double tau_sput = tau_sput_ref * sput_struct;
        if (tau_sput > 0.0 && std::isfinite(tau_sput)) {
          inv_tau_loss += 1.0 / tau_sput;
        }

        double dM = -e_fold_loss_rate(rho_dust, dt, inv_tau_loss);
        if (std::isnan(dM)) {
          std::cout << "dM (species) calculated as NaN, " << dM << std::endl;
          dM = 0.0;
        }
        return dM;
      };

      if (carb_active) {
        destruction_dM_carbon[i] = compute_dM(
            rho_dust_carb_i, shock_factor_carbon);
      }
      if (mg_sil_active) {
        destruction_dM_mg_silicate[i] = compute_dM(
            rho_dust_mg_sil_i, shock_factor_silicate);
      }
      if (fe_sil_active) {
        destruction_dM_fe_silicate[i] = compute_dM(
            rho_dust_fe_sil_i, shock_factor_silicate);
      }
    }
  }
}

// ==========================================
// DUST UPDATE (SPECIES-SPECIFIC: Mg-silicate + Fe-silicate + carbonaceous)
// ==========================================
// Per-channel mass exchange between dust species and their corresponding
// gas-phase reactant pools.  Active when dust_species_track == 1.
//
//   carbon channel:   rho_dust_carbonaceous <-> rho_metal_carbon
//   Mg-sil channel:   rho_dust_mg_silicate <-> {Mg, Si, O} as Mg2SiO4
//   Fe-sil channel:   rho_dust_fe_silicate <-> {Fe, Si, O} as Fe2SiO4
//
// Per-channel pre-cap in absolute mass units replaces the legacy 3-way active[]
// shortfall iteration: growth is capped by the limiting reactant, destruction
// is capped by available dust mass.  No SN injection on this path — host code
// seeds the species directly.
void grackle::impl::dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* growth_dM_mg_silicate,
    const double* growth_dM_fe_silicate, const double* growth_dM_carbon,
    const double* destruction_dM_mg_silicate,
    const double* destruction_dM_fe_silicate,
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
  grackle::impl::View<gr_float***> dust_mg_sil(
      my_fields->dust_density_mg_silicate, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_fe_sil(
      my_fields->dust_density_fe_silicate, my_fields->grid_dimension[0],
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

  const double f_mg_sil_Mg = my_chemistry->dust_mg_silicate_f_Mg;
  const double f_mg_sil_Fe = my_chemistry->dust_mg_silicate_f_Fe;
  const double f_mg_sil_Si = my_chemistry->dust_mg_silicate_f_Si;
  const double f_mg_sil_O  = my_chemistry->dust_mg_silicate_f_O;
  const double f_fe_sil_Mg = my_chemistry->dust_fe_silicate_f_Mg;
  const double f_fe_sil_Fe = my_chemistry->dust_fe_silicate_f_Fe;
  const double f_fe_sil_Si = my_chemistry->dust_fe_silicate_f_Si;
  const double f_fe_sil_O  = my_chemistry->dust_fe_silicate_f_O;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] == MASK_FALSE) continue;

    double rho_gas = d(i, idx_range.j, idx_range.k);
    double rho_dust_mg_sil_i = dust_mg_sil(i, idx_range.j, idx_range.k);
    double rho_dust_fe_sil_i = dust_fe_sil(i, idx_range.j, idx_range.k);
    double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);
    double rho_metal_total = metal(i, idx_range.j, idx_range.k);
    double rho_C  = mC(i, idx_range.j, idx_range.k);
    double rho_O  = mO(i, idx_range.j, idx_range.k);
    double rho_Mg = mMg(i, idx_range.j, idx_range.k);
    double rho_Si = mSi(i, idx_range.j, idx_range.k);
    double rho_Fe = mFe(i, idx_range.j, idx_range.k);
    double dt = dt_value[i];
    if (dt <= 0.0) continue;

    double rho_tracked_before = rho_C + rho_O + rho_Mg + rho_Si + rho_Fe;
    // Untracked gas-phase metals (Al, S, Ne, Na, ...). If host code has seeded
    // per-element fields independently and they exceed rho_metal_total,
    // rho_metal_other would go negative and propagate into rho_metal_new on
    // the rebuild below — clamp at zero to prevent silent mass leakage.
    double rho_metal_other = rho_metal_total - rho_tracked_before;
    if (rho_metal_other < 0.0) rho_metal_other = 0.0;

    // dM > 0 means net growth (gas reactants -> dust);
    // dM < 0 means net destruction (dust -> gas reactants).
    double dM_carb = (growth_dM_carbon[i] + destruction_dM_carbon[i]) * dt;
    double dM_mg_sil = (growth_dM_mg_silicate[i] +
                      destruction_dM_mg_silicate[i]) * dt;
    double dM_fe_sil = (growth_dM_fe_silicate[i] +
                      destruction_dM_fe_silicate[i]) * dt;

    // ---------- Per-channel pre-cap ----------
    // Carbon channel: bounded by rho_C (growth) or rho_dust_carb (destruction).
    // Leave a 10% margin on growth so the gas-C reservoir is not driven to
    // exactly zero in a single subcycle — make_consistent's gas-phase C
    // correction (correctCg = Cg/totalCg) would otherwise zero CI/CII/CO/...
    // and destabilize the metal-chemistry implicit solve that follows.
    if (dM_carb > 0.0) {
      dM_carb = std::min(dM_carb, 0.9 * rho_C);
    } else {
      dM_carb = std::max(dM_carb, -rho_dust_carb_i);
    }

    // Mg- and Fe-silicate channels share Si/O, but Mg and Fe are independent.
    // First cap destruction by each species' own dust mass.
    if (dM_mg_sil < 0.0) {
      dM_mg_sil = std::max(dM_mg_sil, -rho_dust_mg_sil_i);
    }
    if (dM_fe_sil < 0.0) {
      dM_fe_sil = std::max(dM_fe_sil, -rho_dust_fe_sil_i);
    }

    // Then cap positive growth jointly against the gas reservoirs. This
    // preserves Mg-silicate growth when Fe alone limits Fe-silicate, and vice
    // versa, while still preventing double-use of shared Si/O. The 0.9
    // safety margin (same rationale as the carbon channel) leaves headroom
    // in each gas reservoir so make_consistent's per-element corrections
    // don't divide by ~zero in the next subcycle.
    const double sil_grow_safety = 0.9;
    double dM_mg_sil_grow = std::max(dM_mg_sil, 0.0);
    double dM_fe_sil_grow = std::max(dM_fe_sil, 0.0);
    for (int cap_iter = 0; cap_iter < 4; cap_iter++) {
      double scale = 1.0;
      int limiter = -1;  // 0=Mg, 1=Fe, 2=Si, 3=O
      auto consider_element = [&](double rho_X, double c_mg_sil,
                                  double c_fe_sil, int element_id) {
        double need = dM_mg_sil_grow * c_mg_sil +
                      dM_fe_sil_grow * c_fe_sil;
        double budget = sil_grow_safety * rho_X;
        if (need > budget && need > 0.0) {
          double trial = std::max(0.0, budget / need);
          if (trial < scale) {
            scale = trial;
            limiter = element_id;
          }
        }
      };
      consider_element(rho_Mg, f_mg_sil_Mg, f_fe_sil_Mg, 0);
      consider_element(rho_Fe, f_mg_sil_Fe, f_fe_sil_Fe, 1);
      consider_element(rho_Si, f_mg_sil_Si, f_fe_sil_Si, 2);
      consider_element(rho_O,  f_mg_sil_O,  f_fe_sil_O,  3);
      if (limiter < 0 || scale >= 1.0) break;

      if (limiter == 0) {
        if (f_mg_sil_Mg > 0.0) dM_mg_sil_grow *= scale;
        if (f_fe_sil_Mg > 0.0) dM_fe_sil_grow *= scale;
      } else if (limiter == 1) {
        if (f_mg_sil_Fe > 0.0) dM_mg_sil_grow *= scale;
        if (f_fe_sil_Fe > 0.0) dM_fe_sil_grow *= scale;
      } else if (limiter == 2) {
        if (f_mg_sil_Si > 0.0) dM_mg_sil_grow *= scale;
        if (f_fe_sil_Si > 0.0) dM_fe_sil_grow *= scale;
      } else {
        if (f_mg_sil_O > 0.0) dM_mg_sil_grow *= scale;
        if (f_fe_sil_O > 0.0) dM_fe_sil_grow *= scale;
      }
    }
    if (dM_mg_sil > 0.0) {
      dM_mg_sil = dM_mg_sil_grow;
    }
    if (dM_fe_sil > 0.0) {
      dM_fe_sil = dM_fe_sil_grow;
    }

    if (dM_mg_sil < 0.0) {
      dM_mg_sil = std::max(dM_mg_sil, -rho_dust_mg_sil_i);
    } else {
      dM_mg_sil = std::max(dM_mg_sil, 0.0);
    }
    if (dM_fe_sil < 0.0) {
      dM_fe_sil = std::max(dM_fe_sil, -rho_dust_fe_sil_i);
    } else {
      dM_fe_sil = std::max(dM_fe_sil, 0.0);
    }

    // ---------- Apply ----------
    rho_dust_carb_i += dM_carb;
    rho_C           -= dM_carb;

    rho_dust_mg_sil_i += dM_mg_sil;
    rho_dust_fe_sil_i += dM_fe_sil;
    rho_Mg -= dM_mg_sil * f_mg_sil_Mg + dM_fe_sil * f_fe_sil_Mg;
    rho_Fe -= dM_mg_sil * f_mg_sil_Fe + dM_fe_sil * f_fe_sil_Fe;
    rho_Si -= dM_mg_sil * f_mg_sil_Si + dM_fe_sil * f_fe_sil_Si;
    rho_O  -= dM_mg_sil * f_mg_sil_O  + dM_fe_sil * f_fe_sil_O;

    // Floors / NaN guard
    rho_dust_carb_i = std::max(0.0, rho_dust_carb_i);
    rho_dust_mg_sil_i = std::max(0.0, rho_dust_mg_sil_i);
    rho_dust_fe_sil_i = std::max(0.0, rho_dust_fe_sil_i);
    rho_C  = std::max(0.0, rho_C);
    rho_O  = std::max(0.0, rho_O);
    rho_Mg = std::max(0.0, rho_Mg);
    rho_Si = std::max(0.0, rho_Si);
    rho_Fe = std::max(0.0, rho_Fe);

    // Bulk dust = Mg-silicate + Fe-silicate + carbonaceous. The compatibility
    // silicate field is the sum of both silicate endmembers.
    double rho_dust_sil_i = rho_dust_mg_sil_i + rho_dust_fe_sil_i;
    double rho_dust_new = rho_dust_sil_i + rho_dust_carb_i;

    // metal_density_other (= total - C - O - Mg - Si - Fe) is unchanged on
    // the dust loop. Rebuilding total metals from the updated tracked fields
    // preserves consistency even if user-edited silicate fractions do not sum
    // to exactly one.
    double rho_tracked_after = rho_C + rho_O + rho_Mg + rho_Si + rho_Fe;
    double rho_metal_new = rho_metal_other + rho_tracked_after;
    rho_metal_new = std::max(0.0, rho_metal_new);
    double delta_metal_total = rho_metal_new - rho_metal_total;

    // Gas density tracks metals as a subset (dust is not part of `density`).
    rho_gas += delta_metal_total;

    if (std::isnan(rho_dust_new) || std::isnan(rho_metal_new) ||
        std::isnan(rho_gas)) {
      std::cout << "dust_update_species: NaN at cell " << i
                << " dM_carb=" << dM_carb
                << " dM_mg_sil=" << dM_mg_sil
                << " dM_fe_sil=" << dM_fe_sil << std::endl;
      continue;
    }

    if (!dryrun) {
      dust(i, idx_range.j, idx_range.k)      = (gr_float)rho_dust_new;
      dust_sil(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust_sil_i;
      dust_mg_sil(i, idx_range.j, idx_range.k) = (gr_float)rho_dust_mg_sil_i;
      dust_fe_sil(i, idx_range.j, idx_range.k) = (gr_float)rho_dust_fe_sil_i;
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
