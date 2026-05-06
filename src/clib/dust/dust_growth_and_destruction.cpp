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
const double species_t_ref = 50;
const double species_nH_ref = 1.0e3;

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
// DUST GROWTH (SPECIES-SPECIFIC: olivine + pyroxene + carbonaceous)
// ==========================================
// Parallel accretion rates onto pre-existing dust seeds, gated by
// dust_species_track == 1. Carbonaceous growth is rate-limited by gas-phase
// carbon. Olivine growth is rate-limited by the least-available reactant in
// {Mg, Si, Fe, O}; pyroxene growth is limited by {Mg, Si, O} because its Fe
// stoichiometric coefficient is zero.
//
// The species tau_ref values follow Hirashita 2011 MNRAS 416, 1340 section
// 2.6: they are normalized at n_H = 1e3 cm^-3, T = 50 K, S = 0.3,
// a = 0.1 micron, and solar abundance of the relevant key species. This branch
// rescales the paper's S_0.3 factor through dust_growth_sticking_coeff and
// a_0.1 through dust_grainsize / 0.1.
// It therefore computes the density factor from local hydrogen number density,
// not from the bulk SIMBA
// dust_growth_densref parameter. Since this path tracks five silicate
// elements, we apply the same key-species logic separately to olivine and
// pyroxene. This preserves tau_ref for a solar mixture while slowing each
// channel only when a reactant required by that channel is depleted.
void grackle::impl::dust_growth_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* growth_dM_olivine, double* growth_dM_pyroxene,
    double* growth_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_oliv(
      my_fields->dust_density_olivine, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_pyro(
      my_fields->dust_density_pyroxene, my_fields->grid_dimension[0],
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
  // tau_ref values are stored in Gyr (Hirashita 2011 section 2.6).
  double tau_ref_sil = my_chemistry->dust_growth_tauref_silicate * 1e9 *
                       sec_per_year / internalu.tbase1;
  double tau_ref_carb = my_chemistry->dust_growth_tauref_carbon * 1e9 *
                        sec_per_year / internalu.tbase1;
  double sticking_factor =
      (my_chemistry->dust_growth_sticking_coeff > 0.0)
          ? 0.3 / my_chemistry->dust_growth_sticking_coeff
          : huge_value;
  double grain_size_factor =
      (my_chemistry->dust_grainsize > 0.0)
          ? my_chemistry->dust_grainsize / 0.1
          : huge_value;

  double f_ol_Mg = my_chemistry->dust_olivine_f_Mg;
  double f_ol_Fe = my_chemistry->dust_olivine_f_Fe;
  double f_ol_Si = my_chemistry->dust_olivine_f_Si;
  double f_ol_O  = my_chemistry->dust_olivine_f_O;
  double f_py_Mg = my_chemistry->dust_pyroxene_f_Mg;
  double f_py_Fe = my_chemistry->dust_pyroxene_f_Fe;
  double f_py_Si = my_chemistry->dust_pyroxene_f_Si;
  double f_py_O  = my_chemistry->dust_pyroxene_f_O;

  double solar_nonmetal =
      std::max(1.0 - my_chemistry->SolarMetalFractionByMass, tiny_value);
  double solar_H = std::max(my_chemistry->HydrogenFractionByMass *
                            solar_nonmetal, tiny_value);
  double solar_C  = my_chemistry->SolarMetalFractionByMass *
                    solar_frac_C / solar_H;
  double solar_O  = my_chemistry->SolarMetalFractionByMass *
                    solar_frac_O / solar_H;
  double solar_Mg = my_chemistry->SolarMetalFractionByMass *
                    solar_frac_Mg / solar_H;
  double solar_Si = my_chemistry->SolarMetalFractionByMass *
                    solar_frac_Si / solar_H;
  double solar_Fe = my_chemistry->SolarMetalFractionByMass *
                    solar_frac_Fe / solar_H;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    growth_dM_olivine[i] = 0.0;
    growth_dM_pyroxene[i] = 0.0;
    growth_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_oliv_i = dust_oliv(i, idx_range.j, idx_range.k);
      double rho_dust_pyro_i = dust_pyro(i, idx_range.j, idx_range.k);
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

      double accr_struct = (species_nH_ref / nH) *
                           std::pow(species_t_ref / T, 0.5);

      // ---------- Carbonaceous: rate-limited by gas-phase C ----------
      if (rho_C >= metal_gate_threshold * rho_gas &&
          rho_dust_carb_i > 0.0 && solar_C > 0.0 && dt > 0.0) {
        double rho_C_total = rho_C + rho_dust_carb_i;
        double z_C_total_eff = std::max(rho_C_total / rho_H_nuclei,
                                        tiny_value);
        double tau_accr_carb = tau_ref_carb * accr_struct *
                               (solar_C / z_C_total_eff) *
                               grain_size_factor *
                               sticking_factor;
        tau_accr_carb = std::min(std::max(tau_accr_carb, tiny_value),
                                 huge_value);
        double frac_avail = rho_C / (rho_dust_carb_i + rho_C);
        frac_avail = std::clamp(frac_avail, 0.0, 1.0);
        double rate = frac_avail * (rho_dust_carb_i / tau_accr_carb);
        growth_dM_carbon[i] = std::min(rate, rho_C / dt);
      }

      // ---------- Silicates: Choban+2022 key-reactant by species ----------
      // For each required element X, the maximum dust mass that could be
      // assembled from X is rho_X / f_X. The bottleneck element sets the
      // rate, and the solar normalization must use that same element.
      auto silicate_growth_rate = [&](double rho_dust, double f_Mg,
                                      double f_Fe, double f_Si,
                                      double f_O) -> double {
        if (rho_dust <= 0.0 || dt <= 0.0) return 0.0;

        double rho_pool = huge_value;
        double solar_pool = 0.0;
        auto consider_pool = [&](double rho_X, double f_X, double solar_X) {
          if (f_X <= 0.0) return;
          double mass_from_X = rho_X / f_X;
          if (mass_from_X < rho_pool) {
            rho_pool = mass_from_X;
            solar_pool = solar_X / f_X;
          }
        };
        consider_pool(rho_Mg, f_Mg, solar_Mg);
        consider_pool(rho_Fe, f_Fe, solar_Fe);
        consider_pool(rho_Si, f_Si, solar_Si);
        consider_pool(rho_O,  f_O,  solar_O);

        if (rho_pool < metal_gate_threshold * rho_gas || solar_pool <= 0.0) {
          return 0.0;
        }

        double rho_total_pool = rho_pool + rho_dust;
        double z_pool_total_eff = std::max(rho_total_pool / rho_H_nuclei,
                                           tiny_value);
        double tau_accr = tau_ref_sil * accr_struct *
                          (solar_pool / z_pool_total_eff) *
                          grain_size_factor *
                          sticking_factor;
        tau_accr = std::min(std::max(tau_accr, tiny_value), huge_value);
        double frac_avail = rho_pool / (rho_dust + rho_pool);
        frac_avail = std::clamp(frac_avail, 0.0, 1.0);
        double rate = frac_avail * (rho_dust / tau_accr);
        return std::min(rate, rho_pool / dt);
      };

      growth_dM_olivine[i] = silicate_growth_rate(
          rho_dust_oliv_i, f_ol_Mg, f_ol_Fe, f_ol_Si, f_ol_O);
      growth_dM_pyroxene[i] = silicate_growth_rate(
          rho_dust_pyro_i, f_py_Mg, f_py_Fe, f_py_Si, f_py_O);
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
// DUST DESTRUCTION (SPECIES-SPECIFIC: olivine + pyroxene + carbonaceous)
// ==========================================
// SN-shock destruction + thermal sputtering applied independently to each
// dust species, gated by dust_species_track == 1. Carbonaceous (graphite) is
// the shock-vulnerability baseline; olivine and pyroxene use the same
// silicate destruction coefficients. Slavin+2015 Table 2 gives gas-cleared
// masses of 990 Msun for silicates and 600 Msun for carbonaceous grains in
// their standard SNR model, so silicates are destroyed about 1.65x faster
// [REF: Slavin, Dwek, Jones 2015 ApJ 803, 7; Jones+1996 ApJ 469, 740].
// Thermal sputtering uses species-specific tau_ref values
// [REF (silicate): Tsai & Mathews 1995 ApJ 448, 84;
//  REF (carbon, ~2x silicate): Nozawa+2006 ApJ 648, 435]
// with the same Draine & Salpeter 1979 / Tielens+1994 scaling form
// (a/0.1) * (1e-27/rho_gas) * ((2e6/T)^2.5 + 1) used by the bulk path.
// Phase D wires the species outputs into dust_update().
void grackle::impl::dust_destruction_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range,
    const gr_mask_type* itmask, const double* dt_value, const double* t_gas,
    double* destruction_dM_olivine, double* destruction_dM_pyroxene,
    double* destruction_dM_carbon) {

  grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_oliv(
      my_fields->dust_density_olivine, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_pyro(
      my_fields->dust_density_pyroxene, my_fields->grid_dimension[0],
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
  // baseline (1.0); silicate follows the standard Slavin+2015 SNR
  // gas-cleared mass ratio, 990/600 = 1.65.
  const double shock_factor_carbon   = 1.0;
  const double shock_factor_silicate = 1.65;

  // Species-specific thermal sputtering tau_ref (params stored in years)
  double tau_sput_ref_sil = my_chemistry->dust_sputter_tauref_silicate *
                            sec_per_year / internalu.tbase1;
  double tau_sput_ref_carb = my_chemistry->dust_sputter_tauref_carbon *
                             sec_per_year / internalu.tbase1;

  // --- MAIN LOOP ---
  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    destruction_dM_olivine[i] = 0.0;
    destruction_dM_pyroxene[i] = 0.0;
    destruction_dM_carbon[i] = 0.0;

    if (itmask[i] != MASK_FALSE) {
      double rho_gas = d(i, idx_range.j, idx_range.k);
      double rho_dust_oliv_i = dust_oliv(i, idx_range.j, idx_range.k);
      double rho_dust_pyro_i = dust_pyro(i, idx_range.j, idx_range.k);
      double rho_dust_carb_i = dust_carb(i, idx_range.j, idx_range.k);

      bool oliv_active = (rho_dust_oliv_i >= dust_gate_threshold * rho_gas);
      bool pyro_active = (rho_dust_pyro_i >= dust_gate_threshold * rho_gas);
      bool carb_active = (rho_dust_carb_i >= dust_gate_threshold * rho_gas);
      if (!oliv_active && !pyro_active && !carb_active) continue;

      double sne_this = use_sne ? sne(i, idx_range.j, idx_range.k) : 0.0;
      if (rho_gas <= 0.0) continue;

      double temp = std::max(t_gas[i], tiny_value);
      double dt = dt_value[i];
      if (dt <= 0.0) continue;

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
      if (oliv_active) {
        destruction_dM_olivine[i] = compute_dM(
            rho_dust_oliv_i, shock_factor_silicate, tau_sput_ref_sil);
      }
      if (pyro_active) {
        destruction_dM_pyroxene[i] = compute_dM(
            rho_dust_pyro_i, shock_factor_silicate, tau_sput_ref_sil);
      }
    }
  }
}

// ==========================================
// DUST UPDATE (SPECIES-SPECIFIC: olivine + pyroxene + carbonaceous)
// ==========================================
// Per-channel mass exchange between dust species and their corresponding
// gas-phase reactant pools.  Active when dust_species_track == 1.
//
//   carbon channel:   rho_dust_carbonaceous <-> rho_metal_carbon
//   olivine channel:  rho_dust_olivine <-> {Mg, Fe, Si, O} as MgFeSiO4
//   pyroxene channel: rho_dust_pyroxene <-> {Mg, Si, O} as MgSiO3
//
// Per-channel pre-cap in absolute mass units replaces the legacy 3-way active[]
// shortfall iteration: growth is capped by the limiting reactant, destruction
// is capped by available dust mass.  No SN injection on this path — host code
// seeds the species directly.
void grackle::impl::dust_update_species(
    chemistry_data* my_chemistry, grackle_field_data* my_fields,
    InternalGrUnits internalu, IndexRange idx_range, const gr_mask_type* itmask,
    const double* dt_value, const double* growth_dM_olivine,
    const double* growth_dM_pyroxene, const double* growth_dM_carbon,
    const double* destruction_dM_olivine,
    const double* destruction_dM_pyroxene,
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
  grackle::impl::View<gr_float***> dust_oliv(
      my_fields->dust_density_olivine, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> dust_pyro(
      my_fields->dust_density_pyroxene, my_fields->grid_dimension[0],
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

  const double f_ol_Mg = my_chemistry->dust_olivine_f_Mg;
  const double f_ol_Fe = my_chemistry->dust_olivine_f_Fe;
  const double f_ol_Si = my_chemistry->dust_olivine_f_Si;
  const double f_ol_O  = my_chemistry->dust_olivine_f_O;
  const double f_py_Mg = my_chemistry->dust_pyroxene_f_Mg;
  const double f_py_Fe = my_chemistry->dust_pyroxene_f_Fe;
  const double f_py_Si = my_chemistry->dust_pyroxene_f_Si;
  const double f_py_O  = my_chemistry->dust_pyroxene_f_O;

  for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
    if (itmask[i] == MASK_FALSE) continue;

    double rho_gas = d(i, idx_range.j, idx_range.k);
    double rho_dust_oliv_i = dust_oliv(i, idx_range.j, idx_range.k);
    double rho_dust_pyro_i = dust_pyro(i, idx_range.j, idx_range.k);
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
    double rho_metal_other = rho_metal_total - rho_tracked_before;

    // dM > 0 means net growth (gas reactants -> dust);
    // dM < 0 means net destruction (dust -> gas reactants).
    double dM_carb = (growth_dM_carbon[i] + destruction_dM_carbon[i]) * dt;
    double dM_oliv = (growth_dM_olivine[i] +
                      destruction_dM_olivine[i]) * dt;
    double dM_pyro = (growth_dM_pyroxene[i] +
                      destruction_dM_pyroxene[i]) * dt;

    // ---------- Per-channel pre-cap ----------
    // Carbon channel: bounded by rho_C (growth) or rho_dust_carb (destruction).
    if (dM_carb > 0.0) {
      dM_carb = std::min(dM_carb, rho_C);
    } else {
      dM_carb = std::max(dM_carb, -rho_dust_carb_i);
    }

    // Olivine and pyroxene channels share Mg/Si/O but not Fe. First cap
    // destruction by each species' own dust mass.
    if (dM_oliv < 0.0) {
      dM_oliv = std::max(dM_oliv, -rho_dust_oliv_i);
    }
    if (dM_pyro < 0.0) {
      dM_pyro = std::max(dM_pyro, -rho_dust_pyro_i);
    }

    // Then cap positive growth jointly against the shared gas reservoirs.
    // This preserves pyroxene growth when Fe alone limits olivine, while still
    // preventing double-use of Mg/Si/O if both channels grow together.
    double dM_oliv_grow = std::max(dM_oliv, 0.0);
    double dM_pyro_grow = std::max(dM_pyro, 0.0);
    for (int cap_iter = 0; cap_iter < 4; cap_iter++) {
      double scale = 1.0;
      int limiter = -1;  // 0=Mg, 1=Fe, 2=Si, 3=O
      auto consider_element = [&](double rho_X, double c_ol,
                                  double c_py, int element_id) {
        double need = dM_oliv_grow * c_ol + dM_pyro_grow * c_py;
        if (need > rho_X && need > 0.0) {
          double trial = std::max(0.0, rho_X / need);
          if (trial < scale) {
            scale = trial;
            limiter = element_id;
          }
        }
      };
      consider_element(rho_Mg, f_ol_Mg, f_py_Mg, 0);
      consider_element(rho_Fe, f_ol_Fe, f_py_Fe, 1);
      consider_element(rho_Si, f_ol_Si, f_py_Si, 2);
      consider_element(rho_O,  f_ol_O,  f_py_O,  3);
      if (limiter < 0 || scale >= 1.0) break;

      if (limiter == 0) {
        if (f_ol_Mg > 0.0) dM_oliv_grow *= scale;
        if (f_py_Mg > 0.0) dM_pyro_grow *= scale;
      } else if (limiter == 1) {
        if (f_ol_Fe > 0.0) dM_oliv_grow *= scale;
        if (f_py_Fe > 0.0) dM_pyro_grow *= scale;
      } else if (limiter == 2) {
        if (f_ol_Si > 0.0) dM_oliv_grow *= scale;
        if (f_py_Si > 0.0) dM_pyro_grow *= scale;
      } else {
        if (f_ol_O > 0.0) dM_oliv_grow *= scale;
        if (f_py_O > 0.0) dM_pyro_grow *= scale;
      }
    }
    if (dM_oliv > 0.0) {
      dM_oliv = dM_oliv_grow;
    }
    if (dM_pyro > 0.0) {
      dM_pyro = dM_pyro_grow;
    }

    if (dM_oliv < 0.0) {
      dM_oliv = std::max(dM_oliv, -rho_dust_oliv_i);
    } else {
      dM_oliv = std::max(dM_oliv, 0.0);
    }
    if (dM_pyro < 0.0) {
      dM_pyro = std::max(dM_pyro, -rho_dust_pyro_i);
    } else {
      dM_pyro = std::max(dM_pyro, 0.0);
    }

    // ---------- Apply ----------
    rho_dust_carb_i += dM_carb;
    rho_C           -= dM_carb;

    rho_dust_oliv_i += dM_oliv;
    rho_dust_pyro_i += dM_pyro;
    rho_Mg -= dM_oliv * f_ol_Mg + dM_pyro * f_py_Mg;
    rho_Fe -= dM_oliv * f_ol_Fe + dM_pyro * f_py_Fe;
    rho_Si -= dM_oliv * f_ol_Si + dM_pyro * f_py_Si;
    rho_O  -= dM_oliv * f_ol_O  + dM_pyro * f_py_O;

    // Floors / NaN guard
    rho_dust_carb_i = std::max(0.0, rho_dust_carb_i);
    rho_dust_oliv_i = std::max(0.0, rho_dust_oliv_i);
    rho_dust_pyro_i = std::max(0.0, rho_dust_pyro_i);
    rho_C  = std::max(0.0, rho_C);
    rho_O  = std::max(0.0, rho_O);
    rho_Mg = std::max(0.0, rho_Mg);
    rho_Si = std::max(0.0, rho_Si);
    rho_Fe = std::max(0.0, rho_Fe);

    // Bulk dust = olivine + pyroxene + carbonaceous. The compatibility
    // silicate field is olivine + pyroxene.
    double rho_dust_sil_i = rho_dust_oliv_i + rho_dust_pyro_i;
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
                << " dM_oliv=" << dM_oliv
                << " dM_pyro=" << dM_pyro << std::endl;
      continue;
    }

    if (!dryrun) {
      dust(i, idx_range.j, idx_range.k)      = (gr_float)rho_dust_new;
      dust_sil(i, idx_range.j, idx_range.k)  = (gr_float)rho_dust_sil_i;
      dust_oliv(i, idx_range.j, idx_range.k) = (gr_float)rho_dust_oliv_i;
      dust_pyro(i, idx_range.j, idx_range.k) = (gr_float)rho_dust_pyro_i;
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
