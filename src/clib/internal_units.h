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

#ifndef __cplusplus
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

  /// erg*cm^3/s per cooling rate unit (this is defined in our comoving unit
  /// system such that it's the same in the proper & comoving frames. In other
  /// words, it's independent of a).
  ///
  /// This is needed in the vast majority of cases where internal units is
  /// accessed so we precompute it
  double coolunit;

} InternalGrUnits;

/// Return the version of hydrogen mass constant used by the internal units
static inline double internalu_get_mh_(InternalGrUnits internalu) {
  // purely for the sake of consistency, this value of the hydrogen mass is
  // initialized as a floating point literal (with the same precision as
  // gr_float and then it is casted to a double)
  return (double)(mh_grflt);
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
  return internalu.urho*(std::pow(internalu.a_value,3))/mh_local_var;
}

/// calculates coefficients used for computing the Jeans' length
///
/// The main benefit of factoring this out is that we can be sure that we are
/// using a consistent choices in the versions of constants
static inline double internalu_calc_coef_ljeans_(InternalGrUnits internalu,
                                                 double gamma) {
  const double mh_local_var = internalu_get_mh_(internalu);
  return std::sqrt((gamma * pi_fortran_val * kboltz_grflt) /
                   (GravConst_grflt * mh_local_var * internalu.dbase1));
}

/// Construct an instance of InternalGrUnits from the frontend_units
static inline InternalGrUnits new_internalu_(
  const code_units* frontend_units
) {
  // rename the frontend_units object for convenience
  const code_units* my_units = frontend_units;

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

  // initialize output units and copy some stuff from frontend units
  InternalGrUnits internalu;
  internalu.a_value = my_units->a_value;
  internalu.a_units = my_units->a_units;
  internalu.extfields_in_comoving = my_units->comoving_coordinates;

  // this represents the name remapping that would occur when passing
  // co_length_units and co_density_units into a Fortran routine
  internalu.uxyz = co_length_units;
  internalu.urho = co_density_units;

  // store the fundamental units of the unit system
  //
  // TODO: (AFTER FINISIHING TRANSCRIPTION) Consider refactoring this logic for
  // the case with `my_units->comoving_coordinates==0` (in that case, we do an
  // unnecessary round-trip)
  internalu.tbase1   = my_units->time_units;
  internalu.xbase1   = internalu.uxyz/(my_units->a_value*my_units->a_units);    // uxyz is [x]*a      = [x]*[a]*a'
  internalu.dbase1   = internalu.urho*pow((my_units->a_value*my_units->a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '

  // lastly, compute coolunit (make sure we use the correct version of mh)
  const double mh_local_var = internalu_get_mh_(internalu);
  internalu.coolunit = (
    pow(my_units->a_units,5) * pow(internalu.xbase1,2) * pow(mh_local_var,2)
  ) / (pow(internalu.tbase1,3) * internalu.dbase1);

  return internalu;
}

#ifndef __cplusplus
} // extern "C"
#endif /* __cplusplus */

#endif /* INTERNAL_UNITS_HPP */
