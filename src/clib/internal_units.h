// See LICENSE file for license and copyright information

/// @file internal_units.h
/// @brief Declares/Implements the InternalGrUnits struct

/// This file creates the InternalGrUnits type.
/// - it is fully C-compatible, but behaves in a loosely object-oriented way
///   (i.e. it has associated functions that act a little like methods)
/// - the impetus for creating the constructs in this file is to aide with
///   transcription from Fortran to C++
/// - you should think of this type (and its associated methods) as a mechansim
///   for grouping together and encapsulating some common logic that was
///   previously repeated across the Fortran versions of the subroutines.
///
/// As a consequence of this type's role in the transcription process:
/// - the logic was preserved almost **EXACTLY** from the original Fortran files
/// - This is done to maximize consistency and minimize the drift in the results
///   in Grackle's answer tests (modifying the logic could have cascading
///   impacts and plausibly make the results drift outside of tolerances, which
///   is obviously something we want to avoid during transcription)
/// - This means that we may be doing some round-about calculations and using
///   some slightly weird versions of some physical constants
///
/// **AFTER WE COMPLETE TRANSCRIPTION**, we should simplify this logic
#ifndef INTERNAL_UNITS_HPP
#define INTERNAL_UNITS_HPP

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */


#include <math.h>

#include "grackle.h"
#include "grackle_macros.h"
#include "phys_constants.h"

/// Encapsulates Grackle’s Internal Unit System. 
///
/// Grackle internally employs an Enzo-style comoving coordinate system (and 
/// internally tracks various quantities/rates in this unit system).
///
/// This struct is commonly employed by Grackle’s internal functions that do
/// most of the heavy lifting. We coerce the user-specified units (described by
/// the frontend's `code_units` type) near the top of the Grackle’s internal
/// function call-stack, to an equivalent description in terms of Grackle’s
/// comoving unit system (the latter description encapsulated by this object).
///
/// @note
/// Once we finish with transcription, we may wish to:
/// - better unify with the code units frontend
/// - adopt more meaningful variable names (most are preserved from the names
///   used in the routines from before transcription)
/// - consider whether or not we even need to track uxyz and urho (I suspect
///   the only reason that these quantities are passed around is that we use
///   them to recompute the other units)
typedef struct InternalGrUnits{
  /// gives the comoving scale-factor when multiplied by a_units
  double a_value;

  /// this is the factor that is multiplied by `a_value` to get the
  /// cosmological scale factor using the standard definition (a=1 when z=0).
  ///
  /// Recall that Grackle's comoving unit system comes from Enzo. Cosmological
  /// Enzo simulations set `a_units = 1/(1+zinit)`, where `zinit` is the
  /// simulation's initial redshift. This choice meanse that, `a_value=1` for
  /// z=zinit.
  double a_units;

  /// specifies whether user specifies fields in Grackle's comoving coordinates
  /// - while this info has historically been passed deeply into the callstack
  ///   I've confirmed that it's only used for converting the field units
  int extfields_in_comoving;

  /// proper-frame: cm per 1 comoving length unit (depends on scale factor, a)
  double uxyz;
  /// comoving-frame: cm per 1 comoving length unit (independent of a)
  double xbase1;

  /// proper frame: g/cm^3 per 1 comoving density unit (depends on a)
  double urho;
  /// comoving frame: g/cm^3 per 1 comoving density unit (independent of a)
  double dbase1;

  /// seconds per time unit (this is defined in our comoving unit system
  /// such that it's the same in the proper & comoving frames. In other words,
  /// it is independent of a)
  double tbase1;

  /// erg*cm^3/s per cooling rate unit (this is defined in Grackle's comoving
  /// unit system such that it's the same in the proper & comoving frames. In
  /// other words, it's independent of the scale factor, a).
  ///
  /// This is needed in the vast majority of cases where internal units is
  /// accessed so we precompute it
  ///
  /// More Details Adapted from Fortran routines
  /// ==========================================
  /// For completeness, the we include a description of this quantity down
  /// below that was adapted from comments in the original Fortran routines for
  /// computing rate coefficients.
  ///
  /// To understand this description, you need to think about "code units" as
  /// being "dimensionless" (this may not be immediately obvious if you think
  /// about code units in manner consistent with the likes of yt). For a
  /// quantity `q`, we define `q = q~ * [q]`, where
  ///
  /// -  q~ denotes the value in "dimensionless code units"
  /// - [q] denotes the "dimensions" (its a multiplicative factor)
  /// -  q  has cgs units
  ///
  /// The original Fortran routines for computing rate coefficients framed this
  /// quantity as the "dimensions of the cooling coefficients (including
  /// constants)." This description is only meaningful in the context of the
  /// equation we show down.
  /// - Note, this equation includes a factor of rho (total density) in the
  ///   denominator because e is the specific energy (energy/unit mass), not
  ///   energy density (i.e. energy/unit volume).
  ///
  /// Here is the equation:
  ///
  /// ```
  /// δe   = Λ    * n1        * n2        * δt     / dens   / a^3
  /// [e]  = Λ    * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
  /// δe~  = Λ~   * n1~       * n2~       * δt~    / dens~  / a~^3 [~]
  /// ```
  ///
  /// Putting things together, `Λ = [Λ] * Λ~` where,
  /// `[Λ] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]`. We can further
  /// simplify this expression by using both
  ///
  /// - `[e] = ([a]*[x])**2 / [time]**2`,
  /// - `[a] = 1 / (1 + zinit)`, which is defined so that a~ = 1 when z=zinit
  ///
  /// The resulting equation is:
  ///   `[Λ] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)`
  ///
  /// > [!note]
  /// > Some of the energy evolution "coefficients" are only used with a
  /// > single power of number density, n. These do not have the /a^3 factor,
  /// > and they have the units
  /// > ```
  /// > [Λ1] = ([a]**2 * [x]**2 * mh) / [time]**3
  /// >      =  [Λ] * [dens] * [a]**3 / mh
  /// > ```
  /// >
  /// > The difference is accounted for through the `dom` variable in cool.src
  /// >
  /// > Some coefficients are used with 3 powers of n, and they have units
  /// > `[Λ3] = [Λ] * mh / [dens] / [a]**3`
  double coolunit;

  /// constant coefficients (including unit conversions) that must be
  /// multiplied by `internal_energy` (in squared velocity code units) when
  /// calculating temperature (this term neglects contributions from mmw and
  /// adiabatic index). For more information, see the documentation of
  /// `get_temperature_units` in the website docs.
  ///
  /// This isn't really a "unit", but we are including it in this struct
  /// because we historically have called it a unit and it used in the vast
  /// majority of cases where this struct gets used.
  ///
  /// Because of the way that relevant units (namely code velocity) in our
  /// comoving unit system are defined, this is the same in the proper &
  /// comoving frame
  double utem;

  /// Specifies which constant to use for the hydrogen mass.
  /// - A value of 1 indicates that we'll use `mh_grflt` (the constant is
  ///   expressed as a value of type gr_float -- this is consistent with what
  ///   Grackle's Fortran routines have historically done).
  /// - A value of 0 indicates that we'll use `mh` (the constant is always
  ///   expressed as a double -- consistent with what has historically happened
  ///   in Grackle's C functions)
  ///
  /// TODO: Once we finish transcription, we should remove this!
  int use_mh_grflt_;

} InternalGrUnits;

/// Return the version of hydrogen mass constant used by the internal units
static inline double internalu_get_mh_(InternalGrUnits internalu) {
  // purely for the sake of consistency, this value of the hydrogen mass is
  // initialized as a floating point literal (with the same precision as
  // gr_float and then it is casted to a double)
  if (internalu.use_mh_grflt_ == 1) {
    return (double)(mh_grflt);
  } else {
    return mh;
  }
}

/// Return the cm*a_unit/s per 1 velocity unit
///
/// The velocity unit is defined in our comoving unit system such that it's the
/// same in the proper & comoving frames. In other words, it's independent of a
static inline double internalu_get_uvel_(InternalGrUnits internalu) {
  return (internalu.uxyz/internalu.a_value) / internalu.tbase1;
}

/// Return the chunit conversion factor
///
/// This quantity is used internally to help calculate the amount of heating
/// performed contributed by certain chemical reactions.
static inline double internalu_get_chunit_(InternalGrUnits internalu){
  double mh_local_var = internalu_get_mh_(internalu);
  double uvel = internalu_get_uvel_(internalu);

  // before it was moved to this function, there were originally 2 versions of
  // this calculation:
  // -> the older, commented-out version had an extra factor of 2 in the
  //    denominator and had a comment stating "1 eV per H2 formed"
  // -> the newer version, without the factor of 2 (shown below), had a
  //    comment stating "1 eV per REACTION (Feb 2020, Gen Chiaki)"
  return (1.60218e-12)/(uvel*uvel*mh_local_var);
}

/// calculates a standard quantity used throughout the codebase
static inline double internalu_calc_dom_(InternalGrUnits internalu) {
  const double mh_local_var = internalu_get_mh_(internalu);
  return internalu.urho*(pow(internalu.a_value,3))/mh_local_var;
}

/// calculates coefficients used for computing the Jeans' length
///
/// The main benefit of factoring this out is that we can be sure that we are
/// using a consistent choices in the versions of constants
static inline double internalu_calc_coef_ljeans_(InternalGrUnits internalu,
                                                 double gamma) {
  const double mh_local_var = internalu_get_mh_(internalu);
  return sqrt((gamma * pi_fortran_val * kboltz_grflt) /
              (GravConst_grflt * mh_local_var * internalu.dbase1));
}


/// computes the conversion factor that is multiplied by the (non-radiative)
/// rate coefficient to get the value in cgs units. (This unit is defined in
/// Grackle's comoving unit system such that it's the same in the proper &
/// comoving frames. In other words, it's independent of the scale factor, a).
///
/// More Detailed Description (Adapted from Fortran routines)
/// =========================================================
/// For completeness, the we include a description of this quantity down
/// below that was adapted from comments in the original Fortran routines for
/// computing rate coefficients.
///
/// To understand this description, you need to think about "code units" as
/// being "dimensionless" (this may not be immediately obvious if you think
/// about code units in manner consistent with the likes of yt). For a quantity
/// `q`, we define `q = q~ * [q]`, where
///
/// -  q~ denotes the value in "dimensionless code units"
/// - [q] denotes the "dimensions" (its a multiplicative factor)
/// -  q  has cgs units
///
/// The original Fortran routines described this quantity as the "dimensions"
/// of the (non-radiative) rate coefficients or [k]. Note, we have chosen to
/// define the the rate coefficients such that they include the units that
/// convert density to number density. As a result of this choice, the rate
/// equations should look like (in dimensionless units, hence the primes):
///
/// ```
/// d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
/// ```
///
/// where:
///  - k~ is the dimenionless rate coefficient
///  - d0~, d1~, and d2~ are 3 dimensionless species mass densities
///    (importantly, `d = [dens]*d~`)
///  - a~ is the dimensionless expansion coefficient. The standard scale factor
///    (i.e. a=1 for z=0) is given by `a = a~ * [a]`. In slightly more detail,
///    `[a] = 1 / (1 + zinit)`, which is defined so that a~ = 1 when z=zinit.
///
/// To get the definition of [k] let's consider:
///
/// ```
/// rate eqn        :  δn0      = k  * n1        * n2        * δt     / a^3
//  rate eqn units  : [dens]/mh = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
/// rate eqn dimless:  δn0~     = k~ * n1~       * n2~       * δt~    / a~^3
/// ```
///
/// Putting this together we find that `k = [k] * k~`, where
/// `[k] = ( [a]^3 * mh ) / ( [dens] * [time] )  (~)`
///
/// **Reminder:** the number densities in these equations are normalized with
/// [dens] which is not a constant (it has a factor a^3), so the number
/// densities must be converted from comoving to proper.
static inline double internalu_calc_kunit_(InternalGrUnits internalu) {
  double uaye = internalu.a_units;
  double mh_local_var = internalu_get_mh_(internalu);
  return (pow(uaye, 3) * mh) / (internalu.dbase1 * internalu.tbase1);
}

// unused and untested:
//static inline double internalu_calc_kunit_3bdy_(InternalGrUnits internalu) {
//  double uaye = internalu.a_units;
//  double mh_local_uvar = internalu_get_mh_(internalu);
//  double kunit = internalu_calc_kunit_(internalu);
//  return kunit*(pow(uaye, 3) * mh) / internalu.dbase1;
//}

/// Helper function for constructing an instance of InternalGrUnits from the
/// frontend_units
///
/// @note
/// You should call new_internalu_ or new_internalu_legacy_C_convention_
/// instead of calling this function directly
static inline InternalGrUnits new_internalu_helper_(
  const code_units* frontend_units,
  int uses_mh_grflt
) {
  // rename the frontend_units object for convenience
  const code_units* my_units = frontend_units;

  // Part 1: Common logic extracted from solve_chemistry and
  //         calculate_(cooling_time|temperature|dust_temperature) that was
  //         historically used to determine unit-related quantities passed into
  //         the fortran subroutines

  // calculate temperature units (obviously, don't change this until we finish
  // transcription -- under the hood, this calculation implicitly chooses a
  // definition of mh that may differ from the value used elesewhere within the
  // functions in the current file)
  double temperature_units = get_temperature_units(my_units);

  // determine the size of comoving length & comoving density units
  // (that are equivalent to the frontend units), as measured in the proper
  // reference frame
  // - this logic (including variable names) has been preserved from the
  //   standard pattern for handling units pre-transcription
  // - in the future, we can rename some of this
  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    // in this case, the frontend unit-system is exactly the same as grackle's
    // internal comoving unit system
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  } else {
    // in this case, the frontend unit-system uses proper units. Here we need
    // to do some coercion to properly represent the values using Grackle's
    // internal unit system
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  // Part 2: initialize output units and copy some stuff from frontend units
  InternalGrUnits internalu;
  internalu.a_value = my_units->a_value;
  internalu.a_units = my_units->a_units;
  internalu.extfields_in_comoving = my_units->comoving_coordinates;

  // store the temperature_units (this name remapping would historically occur)
  internalu.utem = temperature_units;

  // this represents the name remapping that would occur when passing
  // co_length_units and co_density_units into a Fortran routine
  internalu.uxyz = co_length_units;
  internalu.urho = co_density_units;

  // Part 3: Employ logic from fortran subroutines to compute fundamental unit
  //         basis of the unit-system and store them in internalu. Then compute
  //         the cooling units.

  // TODO: (AFTER FINISIHING TRANSCRIPTION) Consider refactoring this logic for
  // the case with `my_units->comoving_coordinates==0` (in that case, we do an
  // unnecessary round-trip)
  internalu.tbase1   = my_units->time_units;
  internalu.xbase1   = internalu.uxyz/(my_units->a_value*my_units->a_units);    // uxyz is [x]*a      = [x]*[a]*a'
  internalu.dbase1   = internalu.urho*pow((my_units->a_value*my_units->a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '

  // lastly, compute coolunit (make sure we use the correct version of mh)
  internalu.use_mh_grflt_ = uses_mh_grflt;
  const double mh_local_var = internalu_get_mh_(internalu);
  internalu.coolunit = (
    pow(my_units->a_units,5) * pow(internalu.xbase1,2) * pow(mh_local_var,2)
  ) / (pow(internalu.tbase1,3) * internalu.dbase1);

  return internalu;
}

/// Construct an instance of InternalGrUnits from the frontend_units
///
/// This is primarily used in routines transcribed from Fortran
static inline InternalGrUnits new_internalu_(
  const code_units* frontend_units
) {
  return new_internalu_helper_(frontend_units, 1);
}

/// Construct an instance of InternalGrUnits from the frontend_units, while
/// following the conventions for defining the hydrogen that have historically
/// been observed in the C layer
///
/// The goal is to eliminate the distinction between this and new_internalu_
/// after we finish up with transcription
static inline InternalGrUnits new_internalu_legacy_C_(
  const code_units* frontend_units
) {
  return new_internalu_helper_(frontend_units, 0);
}


#ifdef __cplusplus
} // extern "C"
#endif /* __cplusplus */

#endif /* INTERNAL_UNITS_HPP */
