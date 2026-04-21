// See LICENSE file for license and copyright information

/// @file fortran_func_wrappers.hpp
/// @brief Declares wrappers around routines that haven't been transcribed yet
///     from Fortran.
///
/// THIS FILE WILL BE DELETED ONCE WE COMPLETE TRANSCRIPTION
///
/// This file exists to aid the transcription process
/// - This file holds C++ functions that wrap untranscribed Fortran routines.
///   The idea is that these wrapper functions take reduced argument lists
///   (that may includes structs) before converting extended arg lists used
///   within the fortran subroutine.
/// - ideally the reduced argument list used by the wrapper function should
///   roughly approximate the final argument list that the routines will use
///   after transcription is complete.

#ifndef FORTRAN_FUNC_WRAPPERS_HPP
#define FORTRAN_FUNC_WRAPPERS_HPP

#ifndef __cplusplus
#error "This file must be read by a c++ compiler"
#endif

#include "dust/multi_grain_species/calc_grain_size_increment_1d.hpp"

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "internal_types.hpp"
#include "internal_units.h"
#include "LUT.hpp"
#include "opaque_storage.hpp"
#include "utils-cpp.hpp"

#include "step_rate_gauss_seidel.hpp"

// callers of these functions are generally expected to locally shorten the
// namespace name when they call these routines
namespace grackle::impl::fortran_wrapper {

/// wrapper for 1d interpolation
inline double interpolate_1d_g(
  double input1,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 1 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {
  double value;
  FORTRAN_NAME(interpolate_1d_g)(
    &input1, gridDim, gridPar1, &dgridPar1, &dataSize, dataField, &value
  );
  return value;
}

/// wrapper for 2d interpolation
inline double interpolate_2d_g(
  double input1, double input2,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 2 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  gr_i64 dataSize, const double* dataField
) {
  double value;
  FORTRAN_NAME(interpolate_2d_g)(
    &input1, &input2,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2,
    &dataSize, dataField, &value
  );
  return value;
}

/// Helper function used to implement interpolate_3dz_g
inline double interpolate_2df3d_g(
  double input1, double input3,
  const gr_i64* GRIMPL_RESTRICT gridDim, // 3 elements
  const double* GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 index2,
  const double* GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize,
  const double* GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_2df3d_g)(
    &input1, &input3,
    gridDim,
    gridPar1, &dgridPar1, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// Similar to interpolate_3d_g except index2 is calculated ahead of time
/// because it corresponds to redshift, which will not change for the entire
/// grid during the current calculation.
inline double interpolate_3dz_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, gr_i64 index2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
  gr_i64 end_int
) {

  double value;
  FORTRAN_NAME(interpolate_3dz_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &end_int, &value
  );
  return value;

}

/// wraps 3d interpolator functions
inline double interpolate_3d_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_3d_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 4d interpolator functions
inline double interpolate_4d_g(
  double input1, double input2, double input3, double input4,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 4 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_4d_g)(
    &input1, &input2, &input3, &input4,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 5d interpolation routine
inline double interpolate_5d_g(
  double input1, double input2, double input3, double input4, double input5,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 5 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  const double * GRIMPL_RESTRICT gridPar5, double dgridPar5,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_5d_g)(
    &input1, &input2, &input3, &input4, &input5,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4, gridPar5, &dgridPar5,
    &dataSize, dataField, &value
  );
  return value;

}

} // namespace grackle::impl::fortran_wrapper

#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

