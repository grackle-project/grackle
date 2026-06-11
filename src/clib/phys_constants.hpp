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
#include "support/config.hpp"

// reminder macros are totally unaffected by the presence of this namespace
namespace GRIMPL_NAMESPACE_DECL {

/// @brief This namespace specifies physical constants (in CGS units)
///
/// The constants ending with a _grflt underscore exist in order to aide with
/// transcription (the original fortran logic defined macro-constants that
/// expanded to literals of type gr_float). Now that transcription is complete
/// it would be nice to remove all versions of the constants with the _gr_flt
/// suffix (but that will require an update to the gold standard).
namespace constants {

//************************************************
// Physics constants
//************************************************

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

//************************************************
// Astronomical constant
//************************************************

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

//*************************************************
// Miscellaneous values adopted from phys_const.def
//*************************************************

/// @brief Planck's constant [g cm^2 s^-1]
inline constexpr gr_float hplanck_grflt = static_cast<gr_float>(6.6260693e-27);

/// @brief ergs per electronvolt (eV)
inline constexpr double ev2erg = 1.60217653e-12;
inline constexpr gr_float ev2erg_grflt = static_cast<gr_float>(ev2erg);

/// @brief Stefan–Boltzmann constant [g s^-3 K^-4]
inline constexpr gr_float sigma_sb_grflt = static_cast<gr_float>(5.670373e-5);

}  // namespace constants

/// @brief defines constants for converting between density and number density
///
/// End-users of Grackle specify species in terms of density fields related
/// to (but not necessarily equal to) the mass density field of the species.
/// The mass_factor constant is used to convert between that quantity and the
/// number density of the species.
///
/// There is inherently a lot of ambiguity about the precise definition of a
/// mass_factor. This discussion aims formalize a description for an implicit
/// that was hard-coded throughout the codebase.
///
/// What is a mass factor?
/// ======================
/// In more detail:
/// - Consider a species X. If an end user specifies a density ρ, then its
///   number density is ``n = ρ / (mass_factor::X * m_H * grams_per_massUnit)``
/// - ``grams_per_massUnit`` is the conversion factor between grams and the
///   mass code unit.
/// - in this formula, ``m_H`` corresponds to the Hydrogen mass density
///   in units of grams
/// - the use of ``m_H`` in this formula probably seem a little funny
///   - naively, I might choose to replace ``m_H`` with the conversion factor
///     between 1 Da (AKA a dalton or unified atomic mass unit) and grams
///   - however the use of ``m_H`` is currently mandated by our definition
///     of @ref grackle_field_data::e_density
///   - (I suspect it also may implicitly be required by the way that we
///     implement the calculation of various cooling rates)
///
/// mass_factor has historically been an integer value. To be clear, there is
/// no reason requiring this choice (i.e. it would obviously be more correct
/// to use precise values). This choice was presumably made because the
/// mass_factor was hardcoded anywhere.
///
/// Because of mass_factor is an integer, there's an ambiguity about how to
/// compute it. Some ideas include:
/// - the rounded value of `(m_avg / m_H)` where `m_avg` is the average mass
///   of all isotopes represented by a given symbol (weighted by abundance)
/// - the mass number (number of protons and neutrons) of the most common
///   isotope.
/// While this distinction doesn't matter for the primordial nuclides, the
/// distinctions may matter for metals.
///
/// We currently assume that the number of bound electrons in a species has no
/// effect on the mass_factor (we obviously ignore smaller effects like mass
/// reductions from either the binding energy of electrons or the energy in
/// covalent bonds)
///
/// Future Plans
/// ============
/// In the immediate future, the goal is to gradually start migrating to this
/// more standardized approach for specifying mass factor. The easiest way to
/// do that is to add more constants into this namespace.
///
/// In the medium term, the goal is to make these values dynamically queryable.
/// This is important for 2 reasons:
/// - we really should be providing an API to let external codes do this sort
///   of thing (so that they know how exactly to compute density fields)
/// - more importantly, the number of species and reactions has grown a lot.
///   Other than perhaps the simple 6 and 9 primordial species networks (and
///   perhaps the 12 species networks), we probably want to avoid hardcoding
///   logic for every single species. To make calculations more dynamic in
///   these regimes, we probably need to dynamically query the values.
///
/// In the long term, once we are confident, that we converted every last
/// implicit use of mass_factor to use the standardized machinery, we should
/// probably consider adopting a more robust definition.
namespace mass_factor {

/// This is the mass_factor for electrons.
///
/// @note
/// The value is tied to the definition of @ref grackle_field_data::e_density
/// (the fact that its 1 is NOT a typo)
inline constexpr double e = 1.0;

// primordial species:
inline constexpr double H = 1.0;
inline constexpr double He = 4.0;
inline constexpr double H2 = 2.0;
inline constexpr double D = 2.0;
inline constexpr double HD = 3.0;

// obviously, we need to fill more
inline constexpr double C = 12.0;
inline constexpr double O = 16.0;
inline constexpr double Mg = 24.0;
inline constexpr double Al = 27.0;
inline constexpr double Si = 28.0;
inline constexpr double S = 32.0;
inline constexpr double Fe = 56.0;

}  // namespace mass_factor

/** @}*/  // end of group

}  // namespace GRIMPL_NAMESPACE_DECL

#endif  // PHYS_CONSTANTS_HPP
