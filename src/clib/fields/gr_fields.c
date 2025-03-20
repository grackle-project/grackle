// See LICENSE file for license and copyright information

/// @file gr_fields.c
/// @brief Implements the gr_fields type

#include <limits.h> // INT_MAX, SIZE_MAX
#include <stdint.h> // int64_t, uint16_t, uint64_t
#include <stdlib.h> // malloc, free

#include "grackle.h"
#include "fields/experimental_api.h"
#include "grackle_macros.h" // GRACKLE_FREE
#include "fields/internal_field_data.h"
#include "fields/gr_fields_private.h"
#include "FrozenKeyIdxBiMap.h"
#include "status_reporting.h"

// helper function
static int free_helper_(gr_fields* fields, int ignore_field_names)
{
  if (fields->field_names != NULL && !ignore_field_names) {
    drop_FrozenKeyIdxBiMap(fields->field_names);
  }
  GRACKLE_FREE(fields->field_data_alloc_props.ptr);
  GRACKLE_FREE(fields->field_slots);
  free(fields);
  return GR_SUCCESS;
}

int grunstableFields_free(gr_fields* fields) { return free_helper_(fields,0); }

// --------------------------------------------------------------------

/// if cur_name is contained by field_names, update the appropriate entry of
/// field_slots with cur_slot & return 1. Otherwise, return 0
static int try_fill_slot_(const char* cur_name, gr_float** cur_slot,
                          const FrozenKeyIdxBiMap* field_names,
                          gr_float*** field_slots)
{
  uint16_t idx = FrozenKeyIdxBiMap_idx_from_key(field_names, cur_name);
  if (idx == STRU16MAP_INVALID_VAL) { return 0; }
  field_slots[idx] = cur_slot;
  return 1;
}

/// set up the contents of field_slots (which is already allocated)
static int setup_slots_(const FrozenKeyIdxBiMap* field_names,
                        internal_field_data* field_data,
                        gr_float*** field_slots)
{
  int count = 0;
  // here, we use XMacros to call try_fill_slot_ for every member (representing
  // a field) in field_data. We add the return code to count

  #define STRINGIFY_(s) #s
  #define ENTRY(NAME)                                                         \
    count += try_fill_slot_(STRINGIFY_(NAME), &(field_data->NAME),            \
                            field_names, field_slots);
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY
  #undef STRINGIFY

  // check that we filled every slot with a corresponding field_names entry
  const int num_fields = FrozenKeyIdxBiMap_size(field_names);
  if (num_fields != count) {
    // report first field name that could not be filled
    for (int i = 0; i < num_fields; i++) {
      if (field_slots[i] == NULL) {
        const char* name = FrozenKeyIdxBiMap_key_from_idx(field_names, i);
        return GrPrintAndReturnErr("`%s` is not a known field name", name);
      }
    }
    GR_INTERNAL_ERROR("Sanity Check Failure: should be unreachable");
  }
  return GR_SUCCESS;
}

/// perform argument-checking for the gr_fields initializer
static int init_arg_check_(gr_fields** out, FrozenKeyIdxBiMap* field_names,
                           uint64_t flags, uint64_t nelements_per_field)
{

  // we are aborting in cases where there is probably a bug
  if (out == NULL) {
    GR_INTERNAL_ERROR("out should not be a nullptr");
  } else if (field_names == NULL) {
    GR_INTERNAL_ERROR("field_names should not be a nullptr");
  }

  // confirm that flags is valid (for now, there are only 2 allowed values)
  if ((flags != 0) && (flags != GR_FCREATE_MANAGED_FIELDDATA)){
    return GrPrintAndReturnErr("flags is invalid");
  }

  // check the value of nelements_per_field
  const int num_fields = FrozenKeyIdxBiMap_size(field_names);
  if ((flags & GR_FCREATE_MANAGED_FIELDDATA) == 0) {
    if (nelements_per_field != 0) { 
      return GrPrintAndReturnErr(
        "nelements_per_field must be 0 when GR_FCREATE_MANAGED_FIELDDATA is "
        "not specified."
      );
    }
  } else if (nelements_per_field == 0) {
    return GrPrintAndReturnErr(
      "nelements_per_field must not be 0 when GR_FCREATE_MANAGED_FIELDDATA is "
      "specified."
    );
  } else if (nelements_per_field > (SIZE_MAX / num_fields)) {
    return GrPrintAndReturnErr("nelements_per_field is too large");
  }

  return GR_SUCCESS;
}


int grFields_from_field_names_(gr_fields** out, FrozenKeyIdxBiMap* field_names,
                               uint64_t flags, uint64_t nelements_per_field)
{
  int ret = init_arg_check_(out, field_names, flags, nelements_per_field);
  if (ret != GR_SUCCESS) { return ret; }

  struct gr_fields_interface_* fields = (struct gr_fields_interface_*)malloc(
    sizeof(struct gr_fields_interface_)
  );
  // quick initialization to lets us use free_helper_(fields,1) for cleanup if
  // something goes wrong
  fields->field_data_alloc_props.ptr = NULL;
  fields->field_slots = NULL;

  // set up field_names
  fields->field_names = field_names;

  // set up fields->field_data and fields->grid_props (we will adjust the
  // grid_layout at the end if we are managing allocations)
  ret = gr_initialize_field_data(&fields->field_data);
  if (ret != GR_SUCCESS) {
    free_helper_(fields,1);
    return ret;
  }
  fields->field_data.grid_dimension = fields->grid_props.dimension;
  fields->field_data.grid_start = fields->grid_props.start;
  fields->field_data.grid_end = fields->grid_props.end;
  for (int i = 0; i < MAX_GRID_RANK; i++) {
    fields->grid_props.dimension[i] = -1;
    fields->grid_props.start[i] = -1;
    fields->grid_props.end[i] = -1;
  }

  // set up fields->field_slots
  const int num_fields = FrozenKeyIdxBiMap_size(field_names);
  fields->field_slots = (gr_float***)malloc(sizeof(gr_float**) * num_fields);
  ret = setup_slots_(field_names, &fields->field_data, fields->field_slots);
  if (ret != GR_SUCCESS) {
    free_helper_(fields,1);
    return ret;
  }

  // set up the field allocation (or lack thereof)
  fields->field_data_alloc_props.elem_per_field = nelements_per_field;
  if (nelements_per_field == 0) {
    fields->field_data_alloc_props.ptr = NULL; // this is a little redundant
  } else {
    gr_float* alloc_ptr = (gr_float*)malloc(
      sizeof(gr_float) * nelements_per_field * num_fields
    );
    fields->field_data_alloc_props.ptr = alloc_ptr;
    // fill in the field_slots
    for (int i = 0; i < num_fields; i++) {
      *(fields->field_slots[i]) = alloc_ptr + (i * nelements_per_field);
    }
  }

  // as the very last thing we do (once everything else is entirely set up)
  if (nelements_per_field != 0) {
    int64_t extent[1] = {(int64_t)nelements_per_field};
    ret = grunstableFields_set_grid_layout(fields, 1, extent, NULL, NULL);
    if (ret != GR_SUCCESS) {
      free_helper_(fields,1);
      return ret;
    }
  }

  *out = fields;
  return GR_SUCCESS;
}

int grunstableFields_init_from_local(
  gr_fields** ptr,
  const chemistry_data* my_chemistry, const chemistry_data_storage* my_rates,
  uint64_t flags, uint64_t nelements_per_field
)
{
  int num_fields = 0;
  const char** field_list = NULL;
  // we need to insert logic to make field_list into a list of pointers to
  // string literals specifying field names
  // -> maybe we can borrow logic from pygrackle?
  GR_INTERNAL_ERROR("Not Implemented Yet");

  // construct field_names from field_list
  FrozenKeyIdxBiMap* field_names_bimap = NULL;
  int ret = new_FrozenKeyIdxBiMap(&field_names_bimap, field_list, num_fields,
                                  BIMAP_REFS_KEYDATA);
  if (ret != GR_SUCCESS) {
    // do we need to cleanup field_list?
    return ret;
  }

  ret = grFields_from_field_names_(ptr, field_names_bimap, flags,
                                   nelements_per_field);
  if (ret != GR_SUCCESS) {
    drop_FrozenKeyIdxBiMap(field_names_bimap);
    return ret;
  }

  return GR_SUCCESS;
}

// --------------------------------------------------------------------

int grunstableFields_set_grid_layout(gr_fields* fields, int grid_rank,
                                     const int64_t* extent,
                                     const int64_t* start, const int64_t* stop)
{
  // used when start and stop are both NULL
  int64_t fallback_start[3], fallback_stop[3];

  // basic argument checking!
  if ((grid_rank < 1) || (grid_rank > 3)) {
    return GrPrintAndReturnErr("grid_rank must be 1, 2, or 3");
  } else if (extent == NULL) {
    return GrPrintAndReturnErr("extent must not be NULL");
  } else if ((start == NULL) != (stop == NULL)) {
    return GrPrintAndReturnErr(
      "start and stop must both be NULL or must both not be NULL"
    );
  } else if (start == NULL) { // prepare the fallback arrays
    start = fallback_start;
    stop = fallback_stop;
    for (int i = 0; i < grid_rank; i++) {
      fallback_start[i] = 0;
      fallback_stop[i] = extent[i];
    }
  }

  // now we inspect the values of extent, start, and stop
  for (int i = 0; i < grid_rank; i++) {
    if (extent[i] > INT_MAX) {
      return GrPrintAndReturnErr("extent entries can't exceed INT_MAX");
    } else if (stop[i] > extent[i]) {
      return GrPrintAndReturnErr("stop[i] can't exceed extent[i]");
    } else if (stop[i] <= start[i]) {
      return GrPrintAndReturnErr("stop[i] must exceed start[i]");
    } else if (start[i] < 0) {
      return GrPrintAndReturnErr("start entries must be non-negative");
    }
  }

  // perform extra checks when the instance manages allocations
  if (fields->field_data_alloc_props.elem_per_field > 0) {
    uint64_t extent_product = 1;
    const char* err_preamble = "When gr_fields manages memory of slots:";
    for (int i = 0; i < grid_rank; i++) {
      extent_product *= (uint64_t)extent[i];

      // in principle, we could remove these requirements...
      if (start[i] != 0) {
        return GrPrintAndReturnErr("%s start entries must be 0", err_preamble);
      } else if (stop[i] != extent[i]) {
        return GrPrintAndReturnErr(
          "%s stop[i] must be the same as extent[i]", err_preamble
        );
      }
    }

    if (extent_product > fields->field_data_alloc_props.elem_per_field) {
      return GrPrintAndReturnErr(
        "The extent specifies a field with more elements than are allocated"
      );
    }

  }

  // apply the mutations (now that we know that there aren’t any issues)
  fields->field_data.grid_rank=3;
  for (int i=0; i < grid_rank; i++) {
    fields->grid_props.dimension[i] = (int)extent[i];
    fields->grid_props.start[i] = (int)start[i];
    fields->grid_props.end[i] = (int)(stop[i] - 1);
  }

  for (int i = grid_rank; i < 3; i++) {
    fields->grid_props.dimension[i] = 1;
    fields->grid_props.start[i] = 0;
    fields->grid_props.end[i] = 0;
  }

  return GR_SUCCESS;
}

int grunstableFields_set_grid_dx(gr_fields* fields, double dx)
{
  fields->field_data.grid_dx = dx;
  return GR_SUCCESS;
}

// --------------------------------------------------------------------

int grunstableFields_num_fields(const gr_fields* fields, int64_t* n_fields)
{
  *n_fields = (int64_t)FrozenKeyIdxBiMap_size(fields->field_names);
  return GR_SUCCESS;
}


int grunstableFields_name_by_idx(const gr_fields* fields, int64_t idx,
                                 const char** name)
{
  const char* ptr = NULL;
  if (idx <= UINT16_MAX) { // this is for avoiding casting errors
    ptr = FrozenKeyIdxBiMap_key_from_idx(fields->field_names, (uint16_t)idx);
  }
  if (ptr == NULL) { return GrPrintAndReturnErr("idx is too large"); }
  *name = ptr;
  return GR_SUCCESS;
}


int grunstableFields_idx_by_name(const gr_fields* fields, const char* name,
                                 int64_t* idx)
{
  uint16_t tmp_idx = FrozenKeyIdxBiMap_idx_from_key(fields->field_names, name);
  if (tmp_idx == STRU16MAP_INVALID_VAL) {
    *idx = -1;
  } else {
    *idx = (int64_t)tmp_idx;
  }
  return GR_SUCCESS;
}

// --------------------------------------------------------------------

int grunstableFields_store_by_name(gr_fields* fields, const char* name,
                                   gr_float* ptr)
{
  if (fields->field_data_alloc_props.elem_per_field > 0) {
    return GrPrintAndReturnErr(
      "cannot store a pointer in a slot when field_data manages allocations"
    );
  }
  uint16_t idx = FrozenKeyIdxBiMap_idx_from_key(fields->field_names, name);
  if (idx == STRU16MAP_INVALID_VAL) {
    return GrPrintAndReturnErr("`%s` is not a tracked field", name);
  }
  *fields->field_slots[idx] = ptr;
  return GR_SUCCESS;
}

gr_float** grunstableFields_get_slot_by_name(gr_fields* fields,
                                             const char* name)
{
  uint16_t idx = FrozenKeyIdxBiMap_idx_from_key(fields->field_names, name);
  if (idx == STRU16MAP_INVALID_VAL) { return NULL; }
  return fields->field_slots[idx];
}

