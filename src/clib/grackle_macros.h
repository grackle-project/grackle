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

#include <cstdlib>

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

/* Standard definitions (well, fairly standard) */

#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

#ifndef FALSE
#define FALSE     0
#define TRUE      1
#endif

#define FLOAT_UNDEFINED  -99999.0
#define INT_UNDEFINED    -99999

#ifndef tiny
#define tiny 1.0e-20
#endif

#ifndef huge
#define huge 1.0e20
#endif

// the following constants are explicitly defined to always match the values
// historically used in the fortran layer (maybe we can consolidate?)
#ifdef GRACKLE_FLOAT_4
#define tiny_fortran_val 1.0e-20f
#define huge_fortran_val 1.0e20f
#elif defined(GRACKLE_FLOAT_8)
#define tiny_fortran_val 1.0e-20
#define huge_fortran_val 1.0e20
#endif

#define tiny8 1.0e-40
#define huge8 1.0e40

#endif
