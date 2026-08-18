//===----------------------------------------------------------------------===//
//
// See the LICENSE file for license and copyright information
// SPDX-License-Identifier: NCSA AND BSD-3-Clause
//
//===----------------------------------------------------------------------===//
///
/// @file
/// Describes physical constants (in CGS units).
///
/// @par History
/// This file comes from the phys_constants.h header file in Enzo that was
/// originally written by Elizabeth Tasker (May, 2005) and later renamed by
/// Daniel Reynolds
///
//===----------------------------------------------------------------------===//

#ifndef PHYS_CONSTANTS_HPP
#define PHYS_CONSTANTS_HPP

#include "grackle.h"  // gr_float
#include "grackle_float.h"
#include "grackle_macros.h"  // GRFLOAT_C
#include "support/config.hpp"

// reminder macros are totally unaffected by the presence of this namespace
namespace GRIMPL_NAMESPACE_DECL {

/// @defgroup PhysConsts Physical Constants
///
/// This group of entities specify Physical Constants (in CGS units)
///
/// Historically, all constants in this file always expanded to floating point
/// values of type `double`.
/// - In constrast, the "phys_const.def" fortran header always defined
///   macro-constants that expand to floating-point values of type `gr_float`.
/// - To aide with transcribing code from Fortran to C/C++ we have
///   defined versions of most constants in this file that expand to constant
///   of type `gr_float`. These constants have the `_grflt` suffix.
///
/// Ideas for the future:
/// =====================
/// - It would be nice to do away with the alternative versions of these
///   constants
/// - It would be nice to convert all constants from MACROS to c++ constans
///   (reminiscent of the style of the <numbers> header from C++ 20)
///   - it's generally considered better practice to use named constants than
///     macros (it can certainly simplify things with debuggers)
///   - we've already started to do this with the pi constants (this was spurred
///     on by the fact that error arose with g++ if we targetted C++ 20)
/** @{*/

/* Physics constants */

/************************************************/

/* Boltzmann's constant [cm2gs-2K-1] or [ergK-1] */

#define kboltz 1.3806504e-16
#define kboltz_grflt GRFLOAT_C(kboltz)

/* Mass of hydrogen [g] */

#define mh 1.67262171e-24
#define mh_grflt GRFLOAT_C(mh)

namespace constants {

/// @brief mass of an electron [g]
inline constexpr double me = 9.10938215e-28;
inline constexpr gr_float me_grflt = static_cast<gr_float>(me);

/// @brief the mathematical constant π
///
/// If/when we adopt C++ 20, we should consider using the value of
/// std::numbers::pi (we probably want to confirm that the value can't vary
/// between standard library implementations)
inline constexpr double pi = 3.14159265358979323846;

// the following matches the value of `pi_val` from "phys_consts.def"
// TODO: we should unify these definitions with the above value of pi
#ifdef GRACKLE_FLOAT_4
inline constexpr float pi_fortran_val = 3.14159265f;
#else
inline constexpr double pi_fortran_val = 3.141592653589793;
#endif
}  // namespace constants

/************************************************/

/* Astronomical constant */

/************************************************/

/* Speed of light [cms-1] */

#define clight 2.99792458e10
#define clight_grflt GRFLOAT_C(clight)

/* Gravitational constant [cm3g-1s-2]*/

#define GravConst 6.67428e-8
#define GravConst_grflt GRFLOAT_C(GravConst)

/* Solar mass [g] */

#define SolarMass 1.9891e33
#define SolarMass_grflt GRFLOAT_C(SolarMass)

/* Megaparsec [cm] */

#define Mpc 3.0857e24
#define Mpc_grflt GRFLOAT_C(Mpc)

#define kpc 3.0857e21
#define kpc_grflt GRFLOAT_C(kpc)

#define pc 3.0857e18
#define pc_grflt GRFLOAT_C(pc)

/************************************************/

/* Miscellaneous values adopted from phys_const.def */

/************************************************/

#define hplanck_grflt GRFLOAT_C(6.6260693e-27)
#define ev2erg_grflt GRFLOAT_C(1.60217653e-12)
#define sigma_sb_grflt GRFLOAT_C(5.670373e-5)

/** @}*/  // end of group

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // PHYS_CONSTANTS_HPP