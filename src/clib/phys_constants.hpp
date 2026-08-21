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

namespace constants {
/// @brief Boltzmann's constant [cm^2 g s^-2 K^-1] or [erg K^-1]
inline constexpr double kboltz = 1.3806504e-16;
inline constexpr gr_float kboltz_grflt = static_cast<gr_float>(kboltz);

/// @brief mass of hydrogen [g]
inline constexpr double mH = 1.67262171e-24;
inline constexpr gr_float mH_grflt = static_cast<gr_float>(mH);

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

/************************************************/

/* Astronomical constant */

/************************************************/

/// @brief Speed of light [cm s^-1]
inline constexpr double clight = 2.99792458e10;
inline constexpr gr_float clight_grflt = static_cast<gr_float>(clight);

/// @brief Gravitational constant [cm^3 g^-1 s^-2]
inline constexpr double GravConst = 6.67428e-8;
inline constexpr gr_float GravConst_grflt = static_cast<gr_float>(GravConst);

/// @brief Solar mass [g]
inline constexpr double SolarMass = 1.9891e33;
inline constexpr gr_float SolarMass_grflt = static_cast<gr_float>(SolarMass);

/// @brief Megaparsec [cm]
inline constexpr gr_float Mpc = 3.0857e24;
inline constexpr gr_float Mpc_grflt = static_cast<gr_float>(Mpc);

/// @brief kiloparsec [cm]
inline constexpr gr_float kpc = 3.0857e21;
inline constexpr gr_float kpc_grflt = static_cast<gr_float>(kpc);

/// @brief parsec [cm]
inline constexpr gr_float pc = 3.0857e18;
inline constexpr gr_float pc_grflt = static_cast<gr_float>(pc);

/************************************************/

/* Miscellaneous values adopted from phys_const.def */

/************************************************/

/// @brief Planck's constant [g cm^2 s^-1]
inline constexpr gr_float hplanck_grflt = static_cast<gr_float>(6.6260693e-27);

/// @brief ergs per electronvolt (eV)
inline constexpr double ev2erg = 1.60217653e-12;
inline constexpr gr_float ev2erg_grflt = static_cast<gr_float>(ev2erg);

/// @brief Stefan–Boltzmann constant [g s^-3 K^-4]
inline constexpr gr_float sigma_sb_grflt = static_cast<gr_float>(5.670373e-5);

}  // namespace constants

/** @}*/  // end of group

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // PHYS_CONSTANTS_HPP