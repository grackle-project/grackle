// this will be converted to a public header

#ifndef __GRACKLE_EXPERIMENTAL_API_H__
#define __GRACKLE_EXPERIMENTAL_API_H__

#include "grackle.h"

#include <stdint.h> // int64_t, uint64_t, UINT64_C

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// the following is an "opaque type" (that :
// - the definition of the struct will never be in a public header
// - the ONLY way to allocate and delete 
struct gr_fields_interface_;

// in the future, we will define the gr_fields type as:
typedef struct gr_fields_interface_ gr_fields;

// if we want to use bitwise flags,
// - the flage values can't be defined in an enum. There is no way to force an
//   enum to pick a particular represenation in C (before C23) and the
//   implementation is free to represent the values as signed integers. This is
//   problematic since bitwise operations are partially undefined for signed
//   integers.
// - the Fortran bindings would need to define the flags as signed integers
#define GR_FCREATE_MANAGED_FIELDDATA UINT64_C(1)

int grunstableFields_init_from_local(
  gr_fields** ptr,
  const chemistry_data* my_chemistry,
  const chemistry_data_storage* my_rates,
  uint64_t flags,
  uint64_t nelements_per_field
);

int grunstableFields_free(gr_fields*);

int grunstableFields_set_grid_layout(gr_fields* fields, int grid_rank,
                                     const int64_t* extent,
                                     const int64_t* start, const int64_t* stop);

int grunstableFields_set_grid_dx(gr_fields* fields, double dx);

int grunstableFields_num_fields(const gr_fields* fields, int64_t* n_fields);

int grunstableFields_name_by_idx(
  const gr_fields* fields, int64_t idx, const char** name
);

int grunstableFields_idx_by_name(
  const gr_fields* fields, const char* name, int64_t* idx
);

int grunstableFields_store_by_name(
  gr_fields* fields, const char* name, gr_float* ptr
);

gr_float** grunstableFields_get_slot_by_name(
  gr_fields* fields, const char* name
);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* __GRACKLE_EXPERIMENTAL_API_H__ */
