#ifndef __typedefs_h_
#define __typedefs_h_
/***********************************************************************
/
/  MISCELANEOUS TYPEDEFS AND ENUMERATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

/* These are the different types of baryon fields. */

#ifdef SMALL_INTS
typedef gr_int field_type;
typedef gr_int boundary_type;
typedef gr_int gravity_boundary_type;
typedef gr_int interpolation_type;
typedef gr_int hydro_method;
typedef gr_int star_type;
typedef gr_int enum_type;
typedef gr_int staggering;
typedef gr_int fieldtype;
#endif

#ifdef LARGE_INTS
typedef long_int field_type;
typedef long_int boundary_type;
typedef long_int gravity_boundary_type;
typedef long_int interpolation_type;
typedef long_int hydro_method;
typedef long_int star_type;
typedef long_int enum_type;
typedef long_int staggering;
typedef long_int fieldtype;
#endif

#endif
