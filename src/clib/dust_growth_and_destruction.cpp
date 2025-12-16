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
    double* dt_value,
    double* t_gas,
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
    
    double dens_proper = internalu.urho * std::pow(internalu.a_value,3);
    double tau_ref = my_chemistry->dust_growth_tauref * 1e9 * sec_per_year/internalu.tbase1;


    // --- MAIN LOOP ---
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {
        
        double rho_gas      = d(i,idx_range.j,idx_range.k);
        double rho_dust = dust(i,idx_range.j,idx_range.k);
        double rho_metal= metal(i,idx_range.j,idx_range.k);
        double temp     = t_gas[i];
        double dt = dt_value[i];

        double tau_accr0 = tau_ref*(my_chemistry->dust_growth_densref/dens_proper)* std::pow(t_ref/temp,0.5);
        double tau_accr = huge_value;
        tau_accr = std::min(tau_accr0*(my_chemistry->SolarMetalFractionByMass/rho_metal),tau_accr);
        double total_density_init = rho_metal + rho_dust;
        double dM = 0;
        dM = dM + std::min((1 - rho_dust/total_density_init)*(rho_dust/tau_accr)*dt, (total_density_init - rho_dust));
        double dM_tau_accr = dM;

        // recalculate metallicity
        dM = std::max(-1*rho_dust, dM);
        dM = std::min(0.9*rho_metal, dM);
        double dM_conserv = 0.0;
        if (rho_dust >= 0.0) {
            rho_dust = rho_dust + dM;
            rho_metal = rho_metal - dM;
        } else {
            dM_conserv = rho_dust;
            rho_dust = rho_dust - dM_conserv;
            rho_metal = rho_metal + dM_conserv;
        }
        rho_gas = rho_gas + (rho_metal - metal(i,idx_range.j,idx_range.k));
        if (rho_dust < 0) {
            std::exit(21);
        }
        double total_density_final = rho_metal + rho_dust;
        if (std::abs(total_density_final - total_density_init) > 1e-8){
            std::exit(21);
        }
        if (dryrun == false) {
            dust(i,idx_range.j,idx_range.k) = rho_dust;
            metal(i,idx_range.j,idx_range.k) = rho_metal; 
            d(i,idx_range.j,idx_range.k) = rho_gas;
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
    double* dt_value,
    double* t_gas,
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
    grackle::impl::View<gr_float***> sne(
        my_fields->SNe_ThisTimeStep, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    double dens_proper = internalu.urho * std::pow(internalu.a_value,3);

    double Ms100 = 6800.0 * my_chemistry->sne_coeff
                 * (100.0 / my_chemistry->sne_shockspeed)
                 * (100.0 / my_chemistry->sne_shockspeed)
                 * SolarMass / (internalu.urho * std::pow(internalu.uxyz,3));

    // --- MAIN LOOP ---
    for (int i = idx_range.i_start; i <= idx_range.i_end; i++) {

        double rho_gas   = d(i,idx_range.j,idx_range.k);
        double rho_dust  = dust(i,idx_range.j,idx_range.k);
        double rho_metal = metal(i,idx_range.j,idx_range.k);
        double sne_this = sne(i,idx_range.j,idx_range.k);
        double temp      = t_gas[i];
        double dt = dt_value[i];
        double tau_dest = 0;

        double total_density_init = rho_metal + rho_dust;
        double dM = 0;

        // destruction by SN shocks
        if (sne_this <= 0) {
            tau_dest = 1e20;
        } else {
            tau_dest = rho_gas/(Ms100*sne_this*my_chemistry->dust_destruction_eff) * dt;
        }

        // destruction by thermal sputtering
        double tau_sput = 1.7e8 * sec_per_year / internalu.tbase1
                        * (my_chemistry->dust_grainsize/0.1) 
                        * (1.0e-27/(dens_proper * rho_gas)) 
                        * (std::pow((2.0e6/temp),2.5)+1.0);

        double dM_shock = 0.0;
        if (sne_this <= 0) {
            dM_shock = 0.0;
        } else {
            dM_shock = std::min(rho_dust/tau_dest*dt, rho_dust);
        }
        if (dM_shock >= rho_dust) {
            if (dM_shock > rho_dust) {
                std::cout << "WARNING: dM_shock > M_dust SNe shock destruction, " << sne_this << ", " << tau_dest << std::endl;
            }
        } else {
            dM_shock = dM_shock + rho_dust / tau_sput *3.0*dt;
            dM_shock = std::min(dM_shock, rho_dust);
        }
        dM = dM - rho_dust * dM_shock;
        if (std::isnan(dM)) {
            std::cout << "dM calculated as NaN, "<< dM << std::endl;
        }
        // recalculate metallicity
        dM = std::max(-1*rho_dust, dM);
        dM = std::min(0.9*rho_metal, dM);
        double dM_conserv = 0.0;
        if (rho_dust >= 0.0) {
            rho_dust = rho_dust + dM;
            rho_metal = rho_metal - dM;
        } else {
            dM_conserv = rho_dust;
            rho_dust = rho_dust - dM_conserv;
            rho_metal = rho_metal + dM_conserv;
        }
        rho_gas = rho_gas + (rho_metal - metal(i,idx_range.j,idx_range.k));
        if (rho_dust < 0) {
            std::exit(21);
        }
        double total_density_final = rho_metal + rho_dust;
        if (std::abs(total_density_final - total_density_init) > 1e-8){
            std::exit(21);
        }

        if (dryrun == false) {
            dust(i,idx_range.j,idx_range.k) = rho_dust;
            metal(i,idx_range.j,idx_range.k) = rho_metal; 
            d(i,idx_range.j,idx_range.k) = rho_gas;
        }
    }
}
