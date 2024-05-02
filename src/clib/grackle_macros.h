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
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

// define the GR_RESTRICT macro be the restrict keyword introduced in C99
// -> restrict isn't technically part of the C++ standard, but it is commonly
//    provided by c++ compilers (so we try to make use of those alternatives)
#ifdef CONFIG_NO_RESTRICT
  #define GR_RESTRICT /* ... */
#elif !defined(__cplusplus) /* simple case (we are compiling C) */
  #define GR_RESTRICT restrict
#elif __GNUC__
  // C++ compilers other than g++ define this macro. To my knowledge, all of
  // them (e.g. clang++, the new & old intel c++ compilers) define the same
  // the restrict-extension in the same way
  #define GR_RESTRICT __restrict
#else
  #define GR_RESTRICT /* ... */
#endif /* GR_RESTRICT */

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
#define MAX_LINE_LENGTH                   512

#ifndef tiny
#define tiny 1.0e-20
#endif

#ifndef huge
#define huge 1.0e20
#endif

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))
#define COS(X) cos((double) (X))
#define SIN(X) sin((double) (X))

#endif
