/***********************************************************************
/
/ Grackle definitions
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_MACROS_H_
#define __GRACKLE_MACROS_H_

#include "grackle_float.h"

/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

#define GRACKLE_FREE(p)				\
  {						\
    if (p != NULL) {				\
      free(p);					\
      p = NULL;					\
    }						\
  }						\

#ifdef CONFIG_THROW_ABORT
#define GRACKLE_FAIL(A) raise(SIGABRT);
#define GRACKLE_VFAIL(A, ...) raise(SIGABRT);
#else
#define GRACKLE_FAIL(A) throw(GrackleFatalException(A, __FILE__, __LINE__));
#define GRACKLE_VFAIL(format, ...) {snprintf(current_error, 254, format, ##__VA_ARGS__); throw(GrackleFatalException(current_error, __FILE__, __LINE__));}
#endif

/* Fortran name generator (cpp blues) */

#if defined(SUN_OLD)
#define FORTRAN_NAME(NAME) NAME/**/_
#endif

#if defined(IRIS4) || defined(CONVEX) || defined(COMPAQ) || defined(SUN) || defined(LINUX) || defined(IA64) || defined(CRAYX1) || defined(XT3)
#define FORTRAN_NAME(NAME) NAME##_
#endif

#if defined(SPP) || defined(SP2) || defined(BGL)
#define FORTRAN_NAME(NAME) NAME
#endif

#ifdef CONFIG_PFLOAT_16
#define PFORTRAN_NAME(NAME) NAME##_c
#else
#define PFORTRAN_NAME(NAME) FORTRAN_NAME(NAME)
#endif

/* Function macro for gr_float literal */

/// @def GRFLOAT_C(DBL_LITERAL)
/// @brief expands to a a floating point literal having the value specified by
///     it argument and the type `gr_float`. The argument must be a literal
///     with the type `double`
///
/// @par More details
/// This is directly analogous to the `INT32_C(ARG)` or `INTMAX_C(ARG)` macros
/// defined by the standard <stdint.h> header, but it is designed for
/// `gr_float` than rather fixed-size integer types. In more detail:
/// - if `sizeof(gr_float) == sizeof(float)` the macro expands to the input
///   argument with the `f` suffix.
/// - otherwise, the macro expands to the input argument
///
/// @par Concrete Example
/// The snippet, `GRFLOAT_C(1.0)` expands to either `1.0f` or `1.0`.
#ifdef GRACKLE_FLOAT_4
  #define INNER_CONCAT_(A, B) A ## B
  #define GRFLOAT_C(DBL_LITERAL) ( INNER_CONCAT_(DBL_LITERAL, f) )
#elif defined(GRACKLE_FLOAT_8)
  #define GRFLOAT_C(DBL_LITERAL) ( DBL_LITERAL )
#endif

/* HDF5 definitions */

#define HDF5_FILE_I4 H5T_STD_I32BE
#define HDF5_FILE_I8 H5T_STD_I64BE
#define HDF5_FILE_R4 H5T_IEEE_F32BE
#define HDF5_FILE_R8 H5T_IEEE_F64BE
#define HDF5_FILE_B8 H5T_STD_B8BE

#define HDF5_I4  H5T_NATIVE_INT
#define HDF5_I8  H5T_NATIVE_LLONG
#define HDF5_R4  H5T_NATIVE_FLOAT
#define HDF5_R8  H5T_NATIVE_DOUBLE
#define HDF5_R16 H5T_NATIVE_LDOUBLE

/* Precision-dependent definitions */

#ifdef GRACKLE_FLOAT_4
#define FSYM "f"
#define ESYM "e"
#endif

#ifdef GRACKLE_FLOAT_8
#define FSYM "lf"
#define ESYM "le"
#endif

#define GSYM "g"

/* Standard definitions (well, fairly standard) */

#ifndef NULL
#define NULL      0
#endif

// TODO: switch every occurence of FAIL to GR_FAIL
//  -> the use of FAIL conflicts with a macro defined by googletest
//  -> for the moment, we provide a crude hack to work around this
#ifndef SKIP_DEF_FAIL
#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1
#endif /* SKIP_DEF_FAIL */

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999
#define MAX_LINE_LENGTH                   512

#ifndef tiny
#define tiny 1.0e-20
#endif

#ifndef huge
#define huge 1.0e20
#endif

// the following 4 are explicitly defined to always match the values used by
// the fortran layer (in the future, maybe we can consolidate?)
#define tiny_fortran_val GRFLOAT_C(1.0e-20)
#define huge_fortran_val GRFLOAT_C(1.0e20)
#define tiny8 1.0e-40
#define huge8 1.0e40

/* Macro definitions (things C should have) */

#ifndef __cplusplus

// we exclude these from C++ source code because min & max can collide with
// the names of some C++ symbols

// we can replace all occurences of max with fmax in the future
#define max(A,B) ((A) > (B) ? (A) : (B))
// I think we can delete this macro right now (I don't think it's EVER used!)
#define min(A,B) ((A) < (B) ? (A) : (B))

// TODO: remove the following 3 macros (they are NEVER used)
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define COS(X) cos((double) (X))
#define SIN(X) sin((double) (X))
#endif /* end of macro defintions for C code (excluded from C++) */

#define POW(X,Y) pow((double) (X), (double) (Y))

#endif
