#include <cmath>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include "dust_growth_and_destruction.hpp"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

namespace {
    const double k_boltz      = 1.3806504e-16;
    const double m_proton     = 1.67262171e-24;
    const double pi_val       = 3.141592653589793;
    const double sec_per_year = 3.155e7;

    const double tiny_value   = 1.0e-20;
    const double huge_value   = 1.0e+20;

    const double t_ref = 20;

}

// ==========================================
// DUST GROWTH (ACCRETION)
// ==========================================
void grackle::impl::dust_growth(
    chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    InternalGrUnits internalu,
    IndexRange idx_range,
    const gr_mask_type* itmask,
    double* dt_value,
    double* t_gas,
    double* growth_dM)
{
    grackle::impl::View<gr_float***> d(
      my_fields->density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> dust(
      my_fields->dust_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    grackle::impl::View<gr_float***> metal(
      my_fields->metal_density, my_fields->grid_dimension[0],
      my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    double dens_proper = internalu.urho * std::pow(internalu.a_value,3);
    double tau_ref = my_chemistry->dust_growth_tauref * 1e9 * sec_per_year/internalu.tbase1;


    // --- MAIN LOOP ---
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        // Initialize to zero
        growth_dM[i] = 0.0;

        if (itmask[i] != MASK_FALSE) {

            double rho_dust = dust(i,idx_range.j,idx_range.k);
            double rho_metal= metal(i,idx_range.j,idx_range.k);
            double temp     = t_gas[i];
            double dt = dt_value[i];

            double tau_accr0 = tau_ref *
            (my_chemistry->dust_growth_densref / dens_proper) *
            std::pow(t_ref / temp, 0.5);
            double rho_metal_eff = std::max(rho_metal, tiny_value);
            double tau_accr = tau_accr0 * (my_chemistry->SolarMetalFractionByMass / rho_metal_eff);
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
            double dM = std::min(growth_rate, rho_metal/dt);

            // Store the calculated mass change in the output array
            growth_dM[i] = dM;
        }

    }
}

// ==========================================
// 2. DUST DESTRUCTION (SNe + SPUTTERING)
// ==========================================
void grackle::impl::dust_destruction(
    chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    InternalGrUnits internalu,
    IndexRange idx_range,
    const gr_mask_type* itmask,
    double* dt_value,
    double* t_gas,
    double* destruction_dM)
{
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
        my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    double dens_proper = internalu.urho * std::pow(internalu.a_value,3);

    double Ms100 = 6800.0 * my_chemistry->sne_coeff
                 * (100.0 / my_chemistry->sne_shockspeed)
                 * (100.0 / my_chemistry->sne_shockspeed)
                 * SolarMass / (internalu.urho * std::pow(internalu.uxyz,3));

    // --- MAIN LOOP ---
    for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        // Initialize to zero
        destruction_dM[i] = 0.0;

        if (itmask[i] != MASK_FALSE) {

            double rho_gas   = d(i,idx_range.j,idx_range.k);
            double rho_dust  = dust(i,idx_range.j,idx_range.k);
            double sne_this = use_sne ? sne(i,idx_range.j,idx_range.k) : 0.0;
            double temp      = t_gas[i];
            double dt = dt_value[i];
            double tau_dest = 0;

            double dM = 0;
            double dM_shock = 0.0;

            if (use_sne) {
                // destruction by SN shocks
                if (sne_this <= 0) {
                    tau_dest = 1e20;
                    // dM_shock = 0.0;
                } else {
                    tau_dest = rho_gas/(Ms100*sne_this*my_chemistry->dust_destruction_eff) * dt;
                    dM_shock = std::min(rho_dust/tau_dest, rho_dust/dt);
                }
            }

            if (temp >= std::pow(10,5)) {
                // destruction by thermal sputtering
                double tau_sput = 1.7e8 * sec_per_year / internalu.tbase1
                                * (my_chemistry->dust_grainsize/0.1)
                                * (1.0e-27/(dens_proper * rho_gas))
                                * (std::pow((2.0e6/temp),2.5)+1.0);

                if (dM_shock >= rho_dust/dt) {
                    if (dM_shock > rho_dust/dt) {
                        std::cout << "WARNING: dM_shock > M_dust SNe shock destruction, " << sne_this << ", " << tau_dest << std::endl;
                    }
                } else {
                    dM_shock = dM_shock + rho_dust / tau_sput *3.0;
                    dM_shock = std::min(dM_shock, rho_dust/dt);
                }
            }
            //dM = - rho_dust * dM_shock;
            dM = -dM_shock;
            if (std::isnan(dM)) {
                std::cout << "dM calculated as NaN, "<< dM << std::endl;
            }

            // Store the calculated mass change in the output array
            destruction_dM[i] = dM;
        }
    }
}

void grackle::impl::dust_update(
    chemistry_data* my_chemistry,
    grackle_field_data* my_fields,
    InternalGrUnits internalu,
    IndexRange idx_range,
    const gr_mask_type* itmask,
    double* dt_value,
    double* growth_dM,
    double* destruction_dM,
    bool dryrun)
{
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

            double rho_gas   = d(i,idx_range.j,idx_range.k);
            double rho_dust  = dust(i,idx_range.j,idx_range.k);
            double rho_metal = metal(i,idx_range.j,idx_range.k);
            double dt = dt_value[i];

            // Get the total mass change (growth + destruction)

            double dM_total = growth_dM[i] + destruction_dM[i];
            dM_total = dM_total * dt;
            // Apply constraints to dM
            dM_total = std::max(-1*rho_dust, dM_total);
            dM_total = std::min(0.9*rho_metal, dM_total);

            // Apply conservation logic (from original code)
            double dM_conserv = 0.0;
            if (rho_dust >= 0.0) {
                rho_dust = rho_dust + dM_total;
                rho_metal = rho_metal - dM_total;
            } else {
                dM_conserv = rho_dust;
                rho_dust = rho_dust - dM_conserv;
                rho_metal = rho_metal + dM_conserv;
            }

            // Adjust gas density to conserve total mass
            rho_gas = rho_gas + (rho_metal - metal(i,idx_range.j,idx_range.k));

            // Safety checks
            if (rho_dust < 0) {
                fprintf(stderr, "ERROR: Negative dust density at cell %d: rho_dust=%e\n", i, rho_dust);
                std::exit(21);
            }

            // fprintf(stderr,
            //         "internal: dt=%e growth_dM=%.10e destruction_dM=%.10e dM_rate=%.15e gas=%.15e dust=%.15e metal=%.15e consv.=%.15e\n",
            //          dt, growth_dM[i], destruction_dM[i], dM_total, rho_gas, rho_dust, rho_metal, rho_dust+rho_metal);

            // Update the fields
            if (dryrun == false) {
                dust(i,idx_range.j,idx_range.k) = (gr_float)rho_dust;
                metal(i,idx_range.j,idx_range.k) = (gr_float)rho_metal;
                d(i,idx_range.j,idx_range.k) = (gr_float)rho_gas;
            }
        }
    }

}