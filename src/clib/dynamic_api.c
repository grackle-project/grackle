/***********************************************************************
/
/ Defines functions that provide access to fields of the chemistry_data 
/ struct by dynamically specifying strings. If external applications use
/ these functions to access/set fields, rather than directly accessing 
/ the fields of the struct, the applications will be to maintain
/ compatability with multiple versions grackle (as more fields get added 
/ to chemistry_data)
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stddef.h>
#include <string.h>
#include "grackle_chemistry_data.h"

// The following enum is only used internally
enum param_type {INT_PARAM, DOUBLE_PARAM, STRING_PARAM};

// initialize _param_list, which holds an entry for each field in
// chemistry_data. This is only used internally

typedef struct {
  const enum param_type type;
  const size_t offset;
  const char * name;
} param_entry;

static const param_entry _param_list[] =
{
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) \
    { TYPE ## _PARAM, offsetof(chemistry_data, FIELD), #FIELD },
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};

// This is only used internally
static const size_t _n_params = sizeof(_param_list) / sizeof(param_entry);


// define functions for accessing field values
// - local_chemistry_data_access_double
// - local_chemistry_data_access_int
//
// The current implementation has O(N) complexity where N is the number of
// fields of chemistry_data. Since configuration just happens once, this
// probably doesn't need to be very fast...
//
// A faster implementation that does a lot less work at runtime is definitely
// possible (it would just involve more code). 


// retrieves a pointer to the field of my_chemistry that's named ``name``.
//
// This returns a NULL pointer if: my_chemistry is NULL, the field doesn't
// exist, or the field doesn't have the specified type
static void* _get_field_ptr(chemistry_data* my_chemistry, const char* name,
                            enum param_type type)
{
  if (my_chemistry == NULL) { return NULL; }

  for (size_t param_index = 0; param_index < _n_params; param_index++){
    const param_entry* entry = _param_list + param_index;
    if ((strcmp(entry->name, name) == 0) & (entry->type == type)){
      return (void*)( (char*)my_chemistry + entry->offset );
    }
  }
  return NULL;
}

double* local_chemistry_data_access_double(chemistry_data* my_chemistry,
                                           const char* param_name)
{ return (double*)_get_field_ptr(my_chemistry, param_name, DOUBLE_PARAM); }

int* local_chemistry_data_access_int(chemistry_data* my_chemistry,
                                     const char* param_name)
{ return (int*)_get_field_ptr(my_chemistry, param_name, INT_PARAM); }

char** local_chemistry_data_access_string(chemistry_data* my_chemistry,
                                          const char* param_name)
{ return (char**)_get_field_ptr(my_chemistry, param_name, STRING_PARAM); }

// define functions for accessing the names of chemistry_data
// - param_name_double
// - param_name_int
//
// These are primarily needed for testing purposes and can be used for
// serialization. The current implementation is slow. A faster alternative
// would define separate lists for int parameters and double parameters

// returns the name of the ``i``th parameter of the specified type. This returns
// NULL when there are ``i`` or fewer parameters of the specified type
static const char* _param_name(size_t i, enum param_type type)
{
  size_t type_count = 0; // # of parameters of type that have been encountered
  for (size_t param_index = 0; param_index < _n_params; param_index++){
    if (_param_list[param_index].type != type){
      continue;
    } else if (type_count == i){
      return _param_list[param_index].name;
    } else {
      type_count++;
    }
  }
  return NULL;
}

const char* param_name_double(size_t i){ return _param_name(i, DOUBLE_PARAM); }
const char* param_name_int(size_t i){ return _param_name(i, INT_PARAM); }
const char* param_name_string(size_t i){ return _param_name(i, STRING_PARAM); }
