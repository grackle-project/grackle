#include <cmath>
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

// Reference temperature for the bulk (dust_model=1, non-species) accretion
// timescale.
const double bulk_growth_t_ref = 20.0;

// COLIBRE accretion / sputtering reference values
// [REF: Trayford+2026 MNRAS 545, staf2040, section 3.3].
const double colibre_growth_t_ref = 10.0;    // K
const double colibre_growth_nH_ref = 10.0;   // cm^-3
const double colibre_sputter_t_ref = 2.0e6;  // K

const double atomic_C  = 12.01;
const double atomic_O  = 16.00;
const double atomic_Mg = 24.31;
const double atomic_Si = 28.09;
const double atomic_Fe = 55.85;

double limited_expm1(double x) {
  if (!std::isfinite(x) || x <= 0.0) return 0.0;
  return std::expm1(std::min(x, 50.0));
}

// Equivalent per-time rates whose product with dt reproduces one exponential
// e-folding update: rate * dt = rho_dust * (exp(+-dt/tau) - 1).
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

// Gate thresholds for dust evolution: skip cells where dust or metals are a
// negligible fraction of the baryon density. The 1e-9 ratio matches the
// metal-poor threshold used in make_consistent; below it, dust evolution
// would just integrate numerical noise.
const double dust_gate_threshold = 1.0e-9;
const double metal_gate_threshold = 1.0e-9;

// Solar-to-local abundance ratio entering the COLIBRE accretion timescale.
// Degenerate inputs yield an effectively infinite timescale (no growth).
double abundance_ratio(double solar_eps, double local_eps) {
  if (solar_eps <= 0.0 || local_eps <= 0.0) return huge_value;
  return solar_eps / local_eps;
}

// Hydrogen-nuclei mass density census. Counts H from whichever species fields
// are evolved so that the n_H entering the dust rates stays consistent with
// the active network (primordial_chemistry / metal_chemistry); falls back to
// HydrogenFractionByMass * (rho_gas - rho_metal) when no primordial species
// are evolved. Returns 0 for pathological cells (metals >= total density).
struct HNucleiCensus {
  bool use_H, use_H2, use_HD, use_HeH, use_metal_H;
  double hfrac;
  grackle::impl::View<gr_float***> d, metal;
  grackle::impl::View<gr_float***> HI, HII, HM, H2I, H2II, HDI, HDII, HeHII;
  grackle::impl::View<gr_float***> OH, H2O, CH, CH2, OHII, H2OII, H3OII;

  HNucleiCensus(const chemistry_data* my_chemistry,
                const grackle_field_data* my_fields)
      : use_H(my_chemistry->primordial_chemistry > 0),
        use_H2(my_chemistry->primordial_chemistry > 1),
        use_HD(my_chemistry->primordial_chemistry > 2),
        use_HeH(my_chemistry->primordial_chemistry > 3),
        use_metal_H(my_chemistry->metal_chemistry == 1),
        hfrac(my_chemistry->HydrogenFractionByMass),
        d(view_(my_fields, my_fields->density, true)),
        metal(view_(my_fields, my_fields->metal_density, true)),
        HI(view_(my_fields, my_fields->HI_density, use_H)),
        HII(view_(my_fields, my_fields->HII_density, use_H)),
        HM(view_(my_fields, my_fields->HM_density, use_H2)),
        H2I(view_(my_fields, my_fields->H2I_density, use_H2)),
        H2II(view_(my_fields, my_fields->H2II_density, use_H2)),
        HDI(view_(my_fields, my_fields->HDI_density, use_HD)),
        HDII(view_(my_fields, my_fields->HDII_density, use_HeH)),
        HeHII(view_(my_fields, my_fields->HeHII_density, use_HeH)),
        OH(view_(my_fields, my_fields->OH_density, use_metal_H)),
        H2O(view_(my_fields, my_fields->H2O_density, use_metal_H)),
        CH(view_(my_fields, my_fields->CH_density, use_metal_H)),
        CH2(view_(my_fields, my_fields->CH2_density, use_metal_H)),
        OHII(view_(my_fields, my_fields->OHII_density, use_metal_H)),
        H2OII(view_(my_fields, my_fields->H2OII_density, use_metal_H)),
        H3OII(view_(my_fields, my_fields->H3OII_density, use_metal_H)) {}

  double rho_H(int i, int j, int k) const {
    double rho_nonmetal = d(i, j, k) - metal(i, j, k);
    if (rho_nonmetal <= 0.0) return 0.0;
    if (!use_H) return hfrac * rho_nonmetal;

    double rho = HI(i, j, k) + HII(i, j, k);
    if (use_H2) {
      rho += HM(i, j, k) + H2I(i, j, k) + H2II(i, j, k);
    }
    if (use_HD) {
      rho += HDI(i, j, k) / 3.0;
    }
    if (use_HeH) {
      rho += HDII(i, j, k) / 3.0 + HeHII(i, j, k) / 5.0;
    }
    if (use_metal_H) {
      rho += OH(i, j, k) / 17.0 + H2O(i, j, k) * (2.0 / 18.0) +
             CH(i, j, k) / 13.0 + CH2(i, j, k) * (2.0 / 14.0) +
             OHII(i, j, k) / 17.0 + H2OII(i, j, k) * (2.0 / 18.0) +
             H3OII(i, j, k) * (3.0 / 19.0);
    }
    return rho;
  }

 private:
  static grackle::impl::View<gr_float***> view_(
      const grackle_field_data* my_fields, gr_float* ptr, bool active) {
    return grackle::impl::View<gr_float***>(
        active ? ptr : my_fields->density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  }
};

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

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    growth_dM[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_dust = dust(i, idx_range.j, idx_range.k);
      // metal_density holds total gas-phase metals (C, O, and all others)
      double rho_metal = metal(i, idx_range.j, idx_range.k);
      double rho_gas = d(i, idx_range.j, idx_range.k);

      // No metals to accrete onto grains — skip growth
      if (rho_metal < metal_gate_threshold * rho_gas) {
        continue;
      }

      double temp = t_gas[i];
      double dt = dt_value[i];

      // tau_accr = tau_ref * (densref / rho_metal,phys) * sqrt(T_ref / T)
      //          * Z_sun: the density and metallicity scalings combine into
      // the single rho_metal term since rho_metal = Z * rho_gas.
      double tau_accr0 = tau_ref *
                         (my_chemistry->dust_growth_densref / dens_proper) *
                         std::pow(bulk_growth_t_ref / temp, 0.5);
      double rho_metal_eff = std::max(rho_metal, tiny_value);
      double tau_accr =
          tau_accr0 * (my_chemistry->SolarMetalFractionByMass / rho_metal_eff);
      tau_accr = std::clamp(tau_accr, tiny_value, huge_value);

      // Fraction of the combined dust+metal reservoir still in the gas
      // phase; throttles growth as metals deplete.
      double frac_metal_available =
          std::clamp(rho_metal / (rho_dust + rho_metal), 0.0, 1.0);
      double growth_rate = frac_metal_available * (rho_dust / tau_accr);
      growth_dM[i] = std::min(growth_rate, rho_metal / dt);
    }
  }
}

// ==========================================
// DUST GROWTH (SPECIES-SPECIFIC: Mg-silicate + Fe-silicate + carbonaceous)
// ==========================================
// COLIBRE-style accretion [REF: Trayford+2026 MNRAS 545, staf2040, section
// 3.3]. Carbonaceous growth is bottlenecked by gas-phase C. Silicate growth
// is computed once for the total silicate reservoir using the maximum
// depletion factor over O, Si, and Mg+Fe, then split between Mg2SiO4 and
// Fe2SiO4 endmembers according to the local gas-phase Mg/Fe number abundance.
// The dense-gas rate uses a clumped hydrogen density n_H' = C n_H; thermal
// sputtering in dust_destruction_species deliberately uses the unclumped n_H.
// Growth is applied as an exponential e-folding update, returned as an
// equivalent per-time rate for dust_update_species().
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
  HNucleiCensus h_census(my_chemistry, my_fields);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);
  // tau_ref values are stored in Myr (COLIBRE Table 2).
  double tau_ref_sil = my_chemistry->dust_growth_tauref_silicate *
                       sec_per_Myr_local / internalu.tbase1;
  double tau_ref_carb = my_chemistry->dust_growth_tauref_carbon *
                        sec_per_Myr_local / internalu.tbase1;
  // COLIBRE reference values S_ref = 0.3, a_ref = 0.1 micron; non-positive
  // parameters disable growth via an effectively infinite timescale.
  double sticking_factor =
      (my_chemistry->dust_growth_sticking_coeff > 0.0)
          ? 0.3 / my_chemistry->dust_growth_sticking_coeff
          : huge_value;
  double grain_size_factor =
      (my_chemistry->dust_grainsize > 0.0)
          ? my_chemistry->dust_grainsize / 0.1
          : huge_value;

  // Solar abundances (per H mass, then per H nucleus) of the tracked
  // elements, used as the depletion reference in the tau scalings.
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

  // Sub-grid clumping boost: C rises log-linearly in n_H from 1 at nH_min to
  // factor_max at nH_max.
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

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    growth_dM_mg_silicate[i] = 0.0;
    growth_dM_fe_silicate[i] = 0.0;
    growth_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_mg_sil_i = dust_mg_sil(i, idx_range.j, idx_range.k);
      double rho_dust_fe_sil_i = dust_fe_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);

      double rho_C  = mC(i, idx_range.j, idx_range.k);
      double rho_O  = mO(i, idx_range.j, idx_range.k);
      double rho_Mg = mMg(i, idx_range.j, idx_range.k);
      double rho_Si = mSi(i, idx_range.j, idx_range.k);
      double rho_Fe = mFe(i, idx_range.j, idx_range.k);

      double T  = std::max(t_gas[i], tiny_value);
      double dt = dt_value[i];

      double rho_H_nuclei = h_census.rho_H(i, idx_range.j, idx_range.k);
      if (rho_gas <= 0.0 || rho_H_nuclei <= 0.0) {
        continue;
      }
      double nH = rho_H_nuclei * dens_proper / mh;
      if (nH <= 0.0) {
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

      // ---------- Carbonaceous: rate-limited by gas-phase C ----------
      if (rho_C >= metal_gate_threshold * rho_gas &&
          rho_dust_carb_i > 0.0 && dt > 0.0) {
        double tau_accr_carb = tau_ref_carb * accr_struct *
                               abundance_ratio(solar_eps_C, eps_C);
        tau_accr_carb = std::clamp(tau_accr_carb, tiny_value, huge_value);
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
        double tau_accr = tau_ref_sil * accr_struct * sil_ratio;
        tau_accr = std::clamp(tau_accr, tiny_value, huge_value);
        double sil_rate = e_fold_growth_rate(rho_dust_sil_i, dt, tau_accr);

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
          double mg_frac = std::clamp(mg_weight / endmember_weight, 0.0, 1.0);
          growth_dM_mg_silicate[i] = sil_rate * mg_frac;
          growth_dM_fe_silicate[i] = sil_rate * (1.0 - mg_frac);
        }
      }
    }
  }
}

// ==========================================
// DUST DESTRUCTION (SNe + SPUTTERING)
// ==========================================
void grackle::impl::dust_destruction(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, double dt_full, const double* t_gas,
    double* destruction_dM) {
  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
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

  // Gas mass shocked to >= 100 km/s per SN, M_s(100) = 6800 Msun
  // (v_s/100 km/s)^-2 [REF: McKee 1989; used as in Li, Narayanan & Dave
  // 2019], converted to code density * code volume.
  double Ms100 = 6800.0 * my_chemistry->sne_coeff *
                 (100.0 / my_chemistry->sne_shockspeed) *
                 (100.0 / my_chemistry->sne_shockspeed) * SolarMass /
                 (internalu.urho * std::pow(internalu.uxyz, 3));

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
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
      double dM_shock = 0.0;

      if (use_tau_dest) {
        // user-provided destruction timescale
        double tau_dest = tau_dest_field(i, idx_range.j, idx_range.k);
        if (tau_dest <= 0) {
          tau_dest = huge_value;
        }
        dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
      } else if (use_sne && sne_this > 0) {
        // SN shock destruction: sne_rate holds the number of SNe per code
        // volume delivered over the external timestep dt_full, so
        // Ms100 * sne_this / dt_full is the shocked-gas mass rate and
        // tau_dest the time to shock all gas in the cell.
        double tau_dest =
            rho_gas /
            (Ms100 * sne_this * my_chemistry->dust_destruction_eff) * dt_full;
        dM_shock = std::min(rho_dust / tau_dest, rho_dust / dt);
      }

      // Thermal sputtering in hot gas [REF: Tsai & Mathews 1995 ApJ 448,
      // 84]: tau_sp = 0.17 Gyr (a/0.1 um) (rho/1e-27 g cm^-3)^-1
      // [(2e6 K / T)^2.5 + 1]. The factor 3 converts the grain-radius
      // e-folding time to a mass-loss time (m ~ a^3).
      if (temp >= 1.0e5 && dM_shock < rho_dust / dt) {
        double tau_sput = 1.7e8 * sec_per_year / internalu.tbase1 *
                          (my_chemistry->dust_grainsize / 0.1) *
                          (1.0e-27 / (dens_proper * rho_gas)) *
                          (std::pow((2.0e6 / temp), 2.5) + 1.0);
        dM_shock = std::min(dM_shock + 3.0 * rho_dust / tau_sput,
                            rho_dust / dt);
      }

      if (std::isnan(dM_shock)) {
        std::cout << "dust_destruction: dM calculated as NaN" << std::endl;
        dM_shock = 0.0;
      }

      destruction_dM[i] = -dM_shock;
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
// timescale tau_sp = tau_ref (a/0.1) (n_H/cm^-3)^-1 [1 + (T/2e6 K)^-2.5].
// Destruction is applied as exponential decay, returned as an equivalent
// per-time rate for dust_update_species().
void grackle::impl::dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, double dt_full,
    const double* t_gas,
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
  HNucleiCensus h_census(my_chemistry, my_fields);

  double dens_proper = internalu.urho * std::pow(internalu.a_value, 3);

  // Gas mass shocked per SN, M_s(100) = 6800 Msun (v_s/100 km/s)^-2
  // [REF: McKee 1989; used as in Li, Narayanan & Dave 2019].
  double Ms100 = 6800.0 * my_chemistry->sne_coeff *
                 (100.0 / my_chemistry->sne_shockspeed) *
                 (100.0 / my_chemistry->sne_shockspeed) * SolarMass /
                 (internalu.urho * std::pow(internalu.uxyz, 3));

  // Species-specific shock-vulnerability multipliers. Graphite is the
  // baseline (1.0); silicate follows the Slavin+2015 SNR gas-cleared mass
  // ratio, 990/600 = 1.65.
  const double shock_factor_carbon   = 1.0;
  const double shock_factor_silicate = 1.65;

  // COLIBRE/Tsai-Mathews thermal sputtering tau_ref (stored in Myr).
  double tau_sput_ref = my_chemistry->dust_sputter_tauref *
                        sec_per_Myr_local / internalu.tbase1;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    destruction_dM_mg_silicate[i] = 0.0;
    destruction_dM_fe_silicate[i] = 0.0;
    destruction_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_mg_sil_i = dust_mg_sil(i, idx_range.j, idx_range.k);
      double rho_dust_fe_sil_i = dust_fe_sil(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);

      bool mg_sil_active =
          (rho_dust_mg_sil_i >= dust_gate_threshold * rho_gas);
      bool fe_sil_active =
          (rho_dust_fe_sil_i >= dust_gate_threshold * rho_gas);
      bool carb_active = (rho_dust_carb_i >= dust_gate_threshold * rho_gas);
      if (!mg_sil_active && !fe_sil_active && !carb_active) continue;
      if (rho_gas <= 0.0) continue;

      double sne_this = use_sne ? sne(i, idx_range.j, idx_range.k) : 0.0;
      double temp = std::max(t_gas[i], tiny_value);
      double dt = dt_value[i];
      if (dt <= 0.0) continue;

      double rho_H_nuclei = h_census.rho_H(i, idx_range.j, idx_range.k);
      if (rho_H_nuclei <= 0.0) continue;
      double nH = rho_H_nuclei * dens_proper / mh;
      if (nH <= 0.0) continue;

      // Common (species-independent) sputtering structural factor.
      double sput_struct = (my_chemistry->dust_grainsize / 0.1) *
                           (1.0 / nH) *
                           (1.0 + std::pow(temp / colibre_sputter_t_ref,
                                           -2.5));

      // Equivalent per-time destruction rate for one species; its product
      // with dt is an exponential e-folding mass loss.
      auto compute_dM = [&](double rho_dust, double shock_factor) -> double {
        double inv_tau_loss = 0.0;

        if (use_tau_dest) {
          double tau_dest = tau_dest_field(i, idx_range.j, idx_range.k);
          if (tau_dest <= 0) tau_dest = huge_value;
          // Apply species multiplier on top of user-supplied tau_dest:
          // higher shock_factor -> shorter tau -> faster destruction.
          tau_dest = tau_dest / shock_factor;
          if (tau_dest > 0.0 && std::isfinite(tau_dest)) {
            inv_tau_loss += 1.0 / tau_dest;
          }
        } else if (use_sne && sne_this > 0.0) {
          // sne_rate holds the number of SNe per code volume delivered over
          // the external timestep dt_full; Ms100 * sne_this is then the
          // shocked-gas mass per volume in dt_full, and dividing by
          // (rho_gas * dt_full) gives the inverse shock-destruction
          // timescale. The exponential update supplies the dt dependence.
          double inv_tau_shock =
              Ms100 * shock_factor * sne_this *
              my_chemistry->dust_destruction_eff / (rho_gas * dt_full);
          if (inv_tau_shock > 0.0 && std::isfinite(inv_tau_shock)) {
            inv_tau_loss += inv_tau_shock;
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
// Growth is capped by the limiting gas-phase reactant, destruction by the
// available dust mass.  No SN injection on this path — host code seeds the
// species directly.
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
    double rho_metal_other =
        std::max(rho_metal_total - rho_tracked_before, 0.0);

    // dM > 0 means net growth (gas reactants -> dust);
    // dM < 0 means net destruction (dust -> gas reactants).
    double dM_carb = (growth_dM_carbon[i] + destruction_dM_carbon[i]) * dt;
    double dM_mg_sil = (growth_dM_mg_silicate[i] +
                      destruction_dM_mg_silicate[i]) * dt;
    double dM_fe_sil = (growth_dM_fe_silicate[i] +
                      destruction_dM_fe_silicate[i]) * dt;

    // ---------- Per-channel caps ----------
    // Carbon channel: bounded by rho_C (growth) or rho_dust_carb
    // (destruction). Leave a 10% margin on growth so the gas-C reservoir is
    // not driven to exactly zero in a single subcycle — make_consistent's
    // gas-phase C correction (correctCg = Cg/totalCg) would otherwise zero
    // CI/CII/CO/... and destabilize the metal-chemistry implicit solve that
    // follows.
    if (dM_carb > 0.0) {
      dM_carb = std::min(dM_carb, 0.9 * rho_C);
    } else {
      dM_carb = std::max(dM_carb, -rho_dust_carb_i);
    }

    // Mg- and Fe-silicate channels share Si/O, but Mg and Fe are independent.
    // Cap destruction by each species' own dust mass.
    if (dM_mg_sil < 0.0) {
      dM_mg_sil = std::max(dM_mg_sil, -rho_dust_mg_sil_i);
    }
    if (dM_fe_sil < 0.0) {
      dM_fe_sil = std::max(dM_fe_sil, -rho_dust_fe_sil_i);
    }

    // Cap positive growth jointly against the gas reservoirs: each pass
    // rescales only the channels that consume the most-limiting element, so
    // Mg-silicate growth survives when Fe alone limits Fe-silicate (and vice
    // versa) while shared Si/O is never double-used. The 0.9 safety margin
    // (same rationale as the carbon channel) leaves headroom in each gas
    // reservoir so make_consistent's per-element corrections don't divide by
    // ~zero in the next subcycle.
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

    // Rebuild total metals from the updated tracked fields plus the unchanged
    // untracked remainder; this stays consistent even if user-edited silicate
    // fractions do not sum to exactly one.
    double rho_tracked_after = rho_C + rho_O + rho_Mg + rho_Si + rho_Fe;
    double rho_metal_new =
        std::max(rho_metal_other + rho_tracked_after, 0.0);
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

// ==========================================
// DUST UPDATE (BULK)
// ==========================================
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

    // dM_exchange > 0: net growth (metal -> dust); < 0: destruction
    // (dust -> metal). Cap destruction by available dust and growth by 90%
    // of gas-phase metals.
    double dM_exchange = (growth_dM[i] + destruction_dM[i]) * dt;
    dM_exchange = std::max(-rho_dust, dM_exchange);
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
