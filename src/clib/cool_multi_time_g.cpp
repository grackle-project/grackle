//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Implements the cool_multi_time_g function
///
//===----------------------------------------------------------------------===//

// This file was initially generated automatically during conversion of the
// cool_multi_time_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_wrappers.hpp"
#include "index_helper.h"
#include "internal_units.h"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

#include "cool_multi_time_g.h"
#include "scale_fields_g-cpp.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void cool_multi_time_g(
  gr_float* cooltime_data_, int imetal, InternalGrUnits internalu,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, photo_rate_storage my_uvb_rates
)
{
  // shorten `grackle::impl::fortran_wrapper` to `f_wrap` within this function
  namespace f_wrap = ::grackle::impl::fortran_wrapper;

  const grackle_index_helper idx_helper = build_index_helper_(my_fields);

  // Convert densities from comoving to 'proper'
  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,(-3)) );
    grackle::impl::scale_fields_g(imetal, factor, my_chemistry, my_fields);
  }


  OMP_PRAGMA("omp parallel")
  {
    // each OMP thread separately initializes/allocates variables defined in
    // the current scope and then enters the for-loop

    grackle::impl::View<gr_float***> cooltime(cooltime_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

    grackle::impl::GrainSpeciesCollection grain_temperatures =
      grackle::impl::new_GrainSpeciesCollection(my_fields->grid_dimension[0]);

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf =
      grackle::impl::new_Cool1DMultiScratchBuf(my_fields->grid_dimension[0]);
 
    grackle::impl::CoolHeatScratchBuf coolingheating_buf =
      grackle::impl::new_CoolHeatScratchBuf(my_fields->grid_dimension[0]);

    // the following variables aren't embedded because they are structs or are
    // used in a number of different internal routines. Sorting these into
    // additional structs (or leaving them free-standing) will become more
    // obvious as we transcribe more routines.

    std::vector<double> p2d(my_fields->grid_dimension[0]);
    std::vector<double> tgas(my_fields->grid_dimension[0]);
    std::vector<double> mmw(my_fields->grid_dimension[0]);
    std::vector<double> tdust(my_fields->grid_dimension[0]);
    std::vector<double> metallicity(my_fields->grid_dimension[0]);
    std::vector<double> dust2gas(my_fields->grid_dimension[0]);
    std::vector<double> rhoH(my_fields->grid_dimension[0]);
    std::vector<double> edot(my_fields->grid_dimension[0]);

    // Iteration mask for multi_cool
    std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
    std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);

    // The following for-loop is a flattened loop over every k,j combination.
    // OpenMP divides this loop between all threads. Within the loop, we
    // complete calculations for the constructed index-range construct
    // (an index range corresponds to an "i-slice")
    OMP_PRAGMA("omp for")
    for (int t = 0; t < idx_helper.outer_ind_size; t++) {
      const IndexRange idx_range = make_idx_range_(t, &idx_helper);
      const int k = idx_range.k; // use 0-based index
      const int j = idx_range.j; // use 0-based index

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        itmask[i] = MASK_TRUE;
      }

      // Compute the cooling rate
      int dummy_iter_arg=1;

      f_wrap::cool1d_multi_g(
        imetal, idx_range, dummy_iter_arg, edot.data(), tgas.data(),
        mmw.data(), p2d.data(), tdust.data(), metallicity.data(),
        dust2gas.data(), rhoH.data(), itmask.data(), itmask_metal.data(),
        my_chemistry, my_rates, my_fields, my_uvb_rates, internalu,
        grain_temperatures, logTlininterp_buf, cool1dmulti_buf,
        coolingheating_buf
      );

      // Compute the cooling time on the slice
      //   (the gamma used here is the same as used to calculate the pressure
      //    in cool1d_multi_g)

      for (int i = idx_range.i_start; i < idx_range.i_stop; i++) {
        double energy = std::fmax(p2d[i]/(my_chemistry->Gamma-1.),
                                  tiny_fortran_val);
        cooltime(i,j,k) = (gr_float)(energy/edot[i]);
      }

    }

    // cleanup temporaries
    grackle::impl::drop_GrainSpeciesCollection(&grain_temperatures);
    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_Cool1DMultiScratchBuf(&cool1dmulti_buf);
    grackle::impl::drop_CoolHeatScratchBuf(&coolingheating_buf);

  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from 'proper'
  if (internalu.extfields_in_comoving == 1)  {
    gr_float factor = (gr_float)(std::pow(internalu.a_value,3) );
    grackle::impl::scale_fields_g(imetal, factor, my_chemistry, my_fields);
  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
