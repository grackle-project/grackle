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
#include "support/index_helper.hpp"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "internal_types.hpp"
#include "internal_units.hpp"
#include "LUT.hpp"
#include "opaque_storage.hpp"
#include "utils-cpp.hpp"

#include "step_rate_gauss_seidel.hpp"

// callers of these functions are generally expected to locally shorten the
// namespace name when they call these routines
namespace grackle::impl::fortran_wrapper {

/// Performs Gauss-Jordan elimination to solve the specified system of linear
/// equations and inverts the coefficient matrix.
///
/// In more detail, it solves the linear matrix `ax=b`, where `a` is the
/// square coefficient matrix, `b` is the right-hand side vector and `x`
/// is the solution vector
///
/// @param[in]     n The number of linear equations being solved
/// @param[in,out] coef_matrix An n by n column major array that initially
///    specifies the coefficient matrix. This function overwrites this matrix
///    with the inverted matrix.
/// @param[in,out] vec An n element array that initially specifies the
///    right-hand side vector. It's overwritten by the solution vector.
///
/// @retval 0 indicates success
/// @retval 1 indicates that the matrix is singular
///
/// > [!important]
/// > This appears to be taken directly from a routine provided in
/// > "Numerical Recipes" and slightly adapted. The original text allows
/// > `b` to hold multiple right-hand side vectors (and then we compute
/// > solutions for each right-hand side vector at the same time). In
/// > contrast, only allow 1 vector).
/// >
/// > The original Numerical Recipes snippet can be found here
/// > https://phys.uri.edu/nigh/NumRec/bookfpdf/f2-1.pdf
/// >
/// > Code from Numerical Recipes CANNOT be merged into Grackle (the
/// > licensing is fundamentally incompatible)
///
/// @todo
/// We need to replace this for reasons highlighted up above. When we do that,
/// we should account that the only place that calls this routine only needs
/// the solution to the system of equations (i.e. it does not need the
/// inverted matrix). For that reason, we may want to rename this to something
/// a little more generic. We should also consider information highlighted
/// within https://github.com/grackle-project/grackle/issues/255
inline int gaussj_g(int n, double* coef_matrix, double* vector) {
  int ierr;
  FORTRAN_NAME(gaussj_g)(&n, coef_matrix, vector, &ierr);
  return ierr;
}

} // namespace grackle::impl::fortran_wrapper

#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

