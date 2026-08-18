/// @file fortran_func_decls.h
/// @brief Header file defining forward declaring function signatures of
///        Fortran subroutines have not been ported to C++ yet.

#ifndef FORTRAN_FN_DECLARATIONS_HPP
#define FORTRAN_FN_DECLARATIONS_HPP

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "grackle_macros.h" // FORTRAN_NAME
#include "grackle.h"        // gr_float
#include "phys_constants.h" // physical constants
#include <stdint.h>         // int32_t
typedef int32_t gr_mask_type;
typedef long long gr_i64;

#define MASK_TRUE 1
#define MASK_FALSE 0

void FORTRAN_NAME(gaussj_g)(
  int* n, double* a_data_ptr, double* b, int* ierr
);


#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

#endif /* FORTRAN_FN_DECLARATIONS_HPP */
