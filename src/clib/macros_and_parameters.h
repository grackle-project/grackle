#ifndef __macros_and_parameters_h_
#define __macros_and_parameters_h_
/***********************************************************************
/  
/ MACRO DEFINITIONS AND PARAMETERS
/
************************************************************************/

#ifdef CONFIG_THROW_ABORT
#define ENZO_FAIL(A) raise(SIGABRT);
#define ENZO_VFAIL(A, ...) raise(SIGABRT);
#else
#define ENZO_FAIL(A) throw(EnzoFatalException(A, __FILE__, __LINE__));
#define ENZO_VFAIL(format, ...) {snprintf(current_error, 254, format, ##__VA_ARGS__); throw(EnzoFatalException(current_error, __FILE__, __LINE__));}
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

/* Precision-related definitions. */

typedef long long long_int;
typedef long double long_double;
typedef unsigned int unsigned_int;
typedef unsigned long long int unsigned_long_int;

/* Previously in hdf4.h */

typedef float        float32;
typedef double       float64;
typedef long double  float128;

/* Macro definitions for portability */

typedef void           *VOIDP;
typedef int            Eint32;
typedef long long int  Eint64;
typedef float          Eflt32;
typedef double         Eflt64;
typedef long double    Eflt128;
typedef long long int  Elong_int;

typedef int            MPI_Arg;

typedef int            HDF5_hid_t;

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

#if defined(INITS32)
#define inits_type float32
#endif

#if defined(INITS64)
#define inits_type float64
#endif

#ifdef SMALL_INTS
#define Eint int
#define Eunsigned_int unsigned_int
#define ISYM "d"
#define IntDataType MPI_INT
#define HDF5_INT HDF5_I4
#define HDF5_FILE_INT HDF5_FILE_I4
#define nint(A) ( (int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) abs((int) (A))
#define ENPY_INT NPY_INT
#define enpy_int npy_int
#endif

#ifdef LARGE_INTS
#define int long_int // CUDA doesn't like this, and who can blame it?
#define Eint long_int
#define Eunsigned_int unsigned_long_int
#define ISYM "lld"
#define IntDataType MPI_LONG_LONG_INT
#define HDF5_INT HDF5_I8
#define HDF5_FILE_INT HDF5_FILE_I8
#define nint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define nlongint(A) ( (long_int) ((A) + 0.5*sign(A)) )
#define ABS(A) labs((long_int) (A))
#define ENPY_INT NPY_LONG
#define enpy_int npy_long
#endif

#ifdef CONFIG_BFLOAT_4
#define BFLOAT_EPSILON 1e-6f
#define Eflt float
#define FSYM "f"
#define ESYM "e"
#define FloatDataType MPI_FLOAT
#ifdef COMPACT_IO
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R4
#else
#define HDF5_REAL HDF5_R4
#define HDF5_FILE_REAL HDF5_FILE_R8
#endif
#ifdef USE_PYTHON
#define ENPY_BFLOAT NPY_FLOAT
#define enpy_bfloat npy_float
#endif
#endif

#ifdef CONFIG_BFLOAT_8
#define BFLOAT_EPSILON 1e-12f
#define Eflt double
#define FSYM "lf"
#define ESYM "le"
#define FloatDataType MPI_DOUBLE
#define float32 TEMP_HOLD_NAME
#define float double
#define TEMP_HOLD_NAME float32
#define HDF5_REAL HDF5_R8
#define HDF5_FILE_REAL HDF5_FILE_R8
#ifdef USE_PYTHON
#define ENPY_BFLOAT NPY_DOUBLE
#define enpy_bfloat npy_double
#endif
#endif

#ifdef CONFIG_PFLOAT_4
#define PFLOAT_EPSILON 1e-6f
#define FLOAT Eflt32
#define PEXP expf
#define PSYM "f"
#define GSYM "g"
#define GOUTSYM ".8g"
#define MY_MPIFLOAT MPI_FLOAT
#define FLOATDataType MPI_FLOAT
#define HDF5_PREC HDF5_R4
#define HDF5_FILE_PREC HDF5_R4
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_FLOAT
#define enpy_pfloat npy_float
#endif
#endif

#ifdef CONFIG_PFLOAT_8
#define PFLOAT_EPSILON 1e-12f
#define FLOAT double
#define PEXP exp
#define PSYM "lf"
#define GSYM "g"
#define GOUTSYM ".14g"
#define MY_MPIFLOAT MPI_DOUBLE
#define FLOATDataType MPI_DOUBLE
#define HDF5_PREC HDF5_R8
#define HDF5_FILE_PREC HDF5_R8
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_DOUBLE
#define enpy_pfloat npy_double
#endif
#endif

#ifdef CONFIG_PFLOAT_16
#define PFLOAT_EPSILON 1e-16f
#define FLOAT long_double
#define PEXP expl
#define PSYM "Lf"
#define GSYM "g"
#define GOUTSYM ".21Lg"
#define MY_MPIFLOAT MPI_LONG_DOUBLE
#define FLOATDataType MPI_LONG_DOUBLE
#define HDF5_PREC HDF5_R16
#define HDF5_FILE_PREC HDF5_R16
#ifdef USE_PYTHON
#define ENPY_PFLOAT NPY_LONGDOUBLE
#define enpy_pfloat npy_longdouble
#endif
#endif

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

/* Macro definitions (things C should have) */

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#define sign(A)  ((A) >  0  ?  1  : -1 )
#define POW(X,Y) pow((double) (X), (double) (Y))
#define COS(X) cos((double) (X))
#define SIN(X) sin((double) (X))

#endif
