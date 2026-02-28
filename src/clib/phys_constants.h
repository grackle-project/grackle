/***********************************************************************
/
/ Physical constants
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __g_phys_constants_h__
#define __g_phys_constants_h__
/***********************************************************************
/  
/ DEFINITION OF PHYSICAL CONSTANTS
/  
/ written by: Elizabeth Tasker (renamed by Daniel Reynolds)
/ date:       May, 2005 
/
/ Note: CGS units
/
/ Historically, all constants in this file always expanded to floating point
/ values of type `double`.
/ - In constrast, the "phys_const.def" fortran header always defined
/   macro-constants that expand to floating-point values of type `gr_float`.
/ - To aide with transcribing code from Fortran to C/C++ we have
/   defined versions of most constants in this file that expand to constant
/   of type `gr_float`. These constants have the `_grflt` suffix.
/ - In the future, it would be nice to do away with the alternative versions
/   of these constants
/
*********************************************************************/

#include "grackle_float.h"
#include "grackle_macros.h" // GRFLOAT_C

/* Physics constants */

/************************************************/

/* Boltzmann's constant [cm2gs-2K-1] or [ergK-1] */

#define kboltz                          1.3806504e-16
#define kboltz_grflt                    GRFLOAT_C(kboltz)

/* Mass of hydrogen [g] */

#define mh                              1.67262171e-24   
#define mh_grflt                        GRFLOAT_C(mh)

/* Mass of an electron [g] */

#define me                              9.10938215e-28
#define me_grflt                        GRFLOAT_C(me)

/* Pi */

#define pi                              3.14159265358979323846

// the following matches the value of `pi_val` from "phys_consts.def"
#ifdef GRACKLE_FLOAT_4
  #define pi_fortran_val 3.14159265f
#else
  #define pi_fortran_val 3.141592653589793
#endif

/************************************************/

/* Astronomical constant */

/************************************************/

/* Speed of light [cms-1] */ 

#define clight                          2.99792458e10
#define clight_grflt                    GRFLOAT_C(clight)

/* Gravitational constant [cm3g-1s-2]*/

#define GravConst                       6.67428e-8
#define GravConst_grflt                 GRFLOAT_C(GravConst)

/* Solar mass [g] */

#define SolarMass                       1.9891e33
#define SolarMass_grflt                 GRFLOAT_C(SolarMass)

/* Megaparsec [cm] */

#define Mpc                             3.0857e24
#define Mpc_grflt                       GRFLOAT_C(Mpc)

#define kpc                             3.0857e21
#define kpc_grflt                       GRFLOAT_C(kpc)

#define pc                              3.0857e18
#define pc_grflt                        GRFLOAT_C(pc)

/************************************************/

/* Miscellaneous values adopted from phys_const.def */

/************************************************/

#define hplanck_grflt  GRFLOAT_C(6.6260693e-27)
#define ev2erg_grflt   GRFLOAT_C(1.60217653e-12)
#define sigma_sb_grflt GRFLOAT_C(5.670373e-5)

#endif
