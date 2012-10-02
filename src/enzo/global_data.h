/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified:   Robert Harkness
/  date:       February, 2004
/  modified:   Robert Harkness, August 12th 2006
/
/  PURPOSE:
/    This is the global data, which should be held to a minimum.  Any changes
/    in this file require changes in: WriteGlobalData,
/    ReadGlobalData and InitializeNew.  
/    This file is dual-purposed:
/        1) read with    DEFINE_STORAGE defined for the (single) definition
/        2) read without DEFINE_STORAGE defined for external linkage
/
************************************************************************/
#ifndef GLOBAL_DATA_DEFINED__
#define GLOBAL_DATA_DEFINED__

#include <stdio.h>
#ifdef DEFINE_STORAGE
# define EXTERN
#else /* DEFINE_STORAGE */
# define EXTERN extern
#endif

EXTERN char current_error[255];

#endif
