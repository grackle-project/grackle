/***********************************************************************
/
/ Calculate temperature related quantities
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cmath>  // std::fmax
#include <cstdio>
#include <vector>

#include "gas_props.hpp"
#include "grackle.h"
#include "internal_units.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "support/View.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Calculate gas temperature and mmw and pass the values to a callback
///      at each valid location.
///
/// @par History
/// written by: Britton Smith May 2015
/// modified1:  February, 2025 by Matthew Abruzzo; ported to C++
///
/// @param[out] callback Function object invoked at every valid index. This
///     should have the signature
///       `void fn(double Tgas, double mmw, int i, int j, int k)`
///     where `i` is the contiguous index and `k` is the slowest index.
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  cloudy_primordial specifies the cloudy table
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
template<typename Fn>
static int calc_T_related_(Fn callback,
                   const chemistry_data* my_chemistry,
                   cloudy_data cloudy_primordial, grackle_field_data* my_fields,
                   InternalGrUnits internalu)
{
  const int imetal = (my_fields->metal_density != NULL) ? 1 : 0;

  const DensityUnitKind du_kind = (internalu.extfields_in_comoving == 1) ?
    DensityUnitKind::COMOVING : DensityUnitKind::PROPER;

  const IndexHelper idx_helper = build_index_helper_(my_fields);

  OMP_PRAGMA("omp parallel") {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

    // these are used to temporarily hold values from each idx_range
    std::vector<double> tgas(my_fields->grid_dimension[0]);
    std::vector<double> rhoH(my_fields->grid_dimension[0]);
    std::vector<double> mmw(my_fields->grid_dimension[0]);
    std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);

    // The following for-loop is a flattened loop over every k,j combination.
    // OpenMP divides this loop between all threads. Within the loop, we
    // complete calculations for the constructed index-range construct
    // (an index range corresponds to an "i-slice")
    OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      // construct an index-range corresponding to "i-slice"
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);

      // Initialize iteration mask to true for all cells
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;
      }

      // calculate the basic gas properties (tgas, mmw, rhoH)
      basic_gas_props(
          tgas.data(), mmw.data(), rhoH.data(), imetal, itmask.data(),
          my_chemistry, &cloudy_primordial, my_fields, internalu, idx_range,
          du_kind);

      // Record the computed temperature values in the output array
      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        callback(tgas[i], mmw[i], i, idx_range.j, idx_range.k);
      }
    }
  }  // OMP_PRAGMA("omp parallel")

  return GR_SUCCESS;
}

}  // namespace GRIMPL_NAMESPACE_DECL
 
extern "C" int local_calculate_temperature(chemistry_data *my_chemistry,
                                           chemistry_data_storage *my_rates,
                                           code_units *my_units,
                                           grackle_field_data *my_fields,
                                           gr_float *temperature)
{
  if (!my_chemistry->use_grackle) { return GR_SUCCESS; }

  // define a callback function that records T values to temperature
  GRIMPL_NS::FortranView<gr_float***> T_view(
    temperature, my_fields->grid_dimension[0],
    my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  auto callback = [T_view](double T_gas, double mmw, int i, int j, int k) {
    T_view(i, j, k) = static_cast<gr_float>(T_gas);
  };

  return GRIMPL_NS::calc_T_related_(callback, my_chemistry,
                                    my_rates->cloudy_primordial, my_fields,
                                    GRIMPL_NS::new_internalu_(my_units));
}


extern "C" int calculate_temperature(code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *temperature)
{
  if (local_calculate_temperature(grackle_data, &grackle_rates, my_units,
                                  my_fields, temperature) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_calculate_temperature.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}

extern "C" int local_calculate_gamma(chemistry_data *my_chemistry,
                                     chemistry_data_storage *my_rates,
                                     code_units *my_units,
                                     grackle_field_data *my_fields,
                                     gr_float *my_gamma)
{

  if (!my_chemistry->use_grackle) {
    return GR_SUCCESS;
  } else if (my_chemistry->primordial_chemistry <= 1) {
    // If molecular hydrogen is not being used, this is trivial
    // (this should not really be called, but provide it just in case).

    const GRIMPL_NS::IndexHelper ind_helper
        = GRIMPL_NS::build_index_helper_(my_fields);
    for (int outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){
      const GRIMPL_NS::FieldFlatIndexRange range = GRIMPL_NS::inner_flat_range_
          (outer_ind, &ind_helper);
      for (int index = range.start; index <= range.end; index++) {
        my_gamma[index] = my_chemistry->Gamma;
      }
    }
    return GR_SUCCESS;
  } else {  // primordial_chemistry >= 2

    GRIMPL_NS::InternalGrUnits internalu = GRIMPL_NS::new_internalu_(my_units);

    // define a callback function to fill in gamma value at each location
    double temperature_units = internalu.utem;
    GRIMPL_NS::FortranView<const gr_float***> eint(
        my_fields->internal_energy, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    GRIMPL_NS::FortranView<gr_float***> gamma(
        my_gamma, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    auto callback = [=](double T_gas, double mmw, int i, int j, int k) {
      double denom = (temperature_units * mmw) * eint(i, j, k);
      double gamma_minus_one = T_gas / denom;
      gamma(i, j, k) = static_cast<gr_float>(1.0 + gamma_minus_one);
    };
    return GRIMPL_NS::calc_T_related_(callback, my_chemistry,
                                      my_rates->cloudy_primordial, my_fields,
                                      internalu);
  }
}

extern "C" int calculate_gamma(code_units *my_units,
                               grackle_field_data *my_fields,
                               gr_float *my_gamma)
{
  if (local_calculate_gamma(grackle_data, &grackle_rates, my_units,
                            my_fields, my_gamma) != GR_SUCCESS) {
    fprintf(stderr, "Error in local_calculate_gamma.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}

extern "C" int local_calculate_pressure(chemistry_data *my_chemistry,
                                        chemistry_data_storage *my_rates,
                                        code_units *my_units,
                                        grackle_field_data *my_fields,
                                        gr_float *pressure)
{
  // this is only done for historical consistency (I'm not sure we actually
  // want to enforce this minimum)
  constexpr gr_float MIN_PRESSURE = 1.0e-20;

  if (!my_chemistry->use_grackle) {
    return GR_SUCCESS;
  } else if (my_chemistry->primordial_chemistry <= 1) {
    // If molecular hydrogen is not being used, this is trivial
    // (this should not really be called, but provide it just in case).

    const GRIMPL_NS::IndexHelper ind_helper
        = GRIMPL_NS::build_index_helper_(my_fields);

    double gm1 = my_chemistry->Gamma - 1.0;
    const gr_float* rho = my_fields->density;
    const gr_float* eint = my_fields->internal_energy;

    // parallelize the k and j loops with OpenMP
    // (these loops are flattened them for better parallelism)
# ifdef _OPENMP
# pragma omp parallel for schedule( runtime )
# endif
    for (int outer_ind = 0; outer_ind < ind_helper.outer_ind_size; outer_ind++){
      GRIMPL_NS::FieldFlatIndexRange range = GRIMPL_NS::inner_flat_range_(
          outer_ind, &ind_helper);

      for (int index = range.start; index <= range.end; index++) {
        double p = gm1 * (rho[index] * eint[index]);
        pressure[index] = std::fmax(static_cast<gr_float>(p), MIN_PRESSURE);
      }
    }
    return GR_SUCCESS;
  } else { // primordial_chemistry >= 2
    GRIMPL_NS::InternalGrUnits internalu = GRIMPL_NS::new_internalu_(my_units);

    // define a callback function to fill in pressure value at each location
    double mH_div_kboltz = internalu.utem;  // <- in code units
    GRIMPL_NS::FortranView<const gr_float***> rho(
        my_fields->density, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
    GRIMPL_NS::FortranView<gr_float***> pressure_view(
        pressure, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    auto callback = [=](double T_gas, double mmw, int i, int j, int k) {
      // here's the basic algebra
      //      P    = ρ * kboltz * T / (mH * μ)
      //    P / ρ  = kboltz * T / (mH * μ)
      //    P / ρ  = T / ((mH / kboltz) * μ)
      double pressure_div_rho = T_gas / (mH_div_kboltz * mmw);
      gr_float p = static_cast<gr_float>(pressure_div_rho * rho(i,j,k));
      pressure_view(i, j, k) = std::fmax(p, MIN_PRESSURE);
    };
    return GRIMPL_NS::calc_T_related_(callback, my_chemistry,
                                      my_rates->cloudy_primordial, my_fields,
                                      internalu);
  }
}

extern "C" int calculate_pressure(code_units *my_units,
                                  grackle_field_data *my_fields,
                                  gr_float *pressure)
{
  if (local_calculate_pressure(grackle_data, &grackle_rates, my_units,
                               my_fields, pressure) != GR_SUCCESS) {
    std::fprintf(stderr, "Error in local_calculate_pressure.\n");
    return GR_FAIL;
  }
  return GR_SUCCESS;
}
 