/***********************************************************************
/
/ Calculate temperature field
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cstdio>
#include <vector>

#include "gas_props.hpp"
#include "grackle.h"
#include "internal_units.hpp"
#include "scale_fields.hpp"
#include "support/config.hpp"
#include "support/index_helper.hpp"
#include "utils-cpp.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GRIMPL_NAMESPACE_DECL {

/// @brief Fill temperature on a 3D grid.
///
/// @par History
/// written by: Britton Smith May 2015
/// modified1:  February, 2025 by Matthew Abruzzo; ported to C++
///
/// @param[out] temperature_data Array where computed values are written
/// @param[in]  my_chemistry specifies various properties
/// @param[in]  cloudy_primordia specifies the cloudy table
/// @param[in]  my_fields specifies all of the field data
/// @param[in]  internalu Specifies unit information
static int calc_T_(gr_float* temperature_data_,
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

    FortranView<gr_float***> temperature(
        temperature_data_, my_fields->grid_dimension[0],
        my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

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
        temperature(i, idx_range.j, idx_range.k) = tgas[i];
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
  return GRIMPL_NS::calc_T_(temperature, my_chemistry,
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
