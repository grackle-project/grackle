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

inline int my_dummy_fn() { return 3; }

} // namespace grackle::impl::fortran_wrapper

#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

