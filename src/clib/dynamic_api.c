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

// initialize _int_param_l, _double_param_l & _string_param_l. They are sized
// lists that respectively hold entries for each int, double, and string field
// in chemistry_data. These are only used internally
typedef struct { const size_t offset; const char * name; } param_entry;
typedef struct { const size_t len; const param_entry * entries; } param_list;

#define INIT_PAR_LIST(ENTRIES) {sizeof(ENTRIES) / sizeof(param_entry), ENTRIES}

#define _LIST_INT_INT(FIELD) {offsetof(chemistry_data, FIELD), #FIELD},
#define _LIST_INT_DOUBLE(FIELD) /* ... */
#define _LIST_INT_STRING(FIELD) /* ... */

#define _LIST_DOUBLE_INT(FIELD) /* ... */
#define _LIST_DOUBLE_DOUBLE(FIELD) {offsetof(chemistry_data, FIELD), #FIELD},
#define _LIST_DOUBLE_STRING(FIELD) /* ... */

#define _LIST_STRING_INT(FIELD) /* ... */
#define _LIST_STRING_DOUBLE(FIELD) /* ... */
#define _LIST_STRING_STRING(FIELD) {offsetof(chemistry_data, FIELD), #FIELD},

static const param_entry _int_param_entries[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) _LIST_INT_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
static const param_list _int_param_l = INIT_PAR_LIST(_int_param_entries);

static const param_entry _double_param_entries[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) _LIST_DOUBLE_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
static const param_list _double_param_l = INIT_PAR_LIST(_double_param_entries);

static const param_entry _string_param_entries[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) _LIST_STRING_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
static const param_list _string_param_l = INIT_PAR_LIST(_string_param_entries);


// define functions for accessing field values
// - local_chemistry_data_access_int
// - local_chemistry_data_access_double
// - local_chemistry_data_access_string
//
// The current implementation has O(N) complexity where N is the number of
// fields of a given type in chemistry_data. Since configuration just happens
// once, this probably doesn't need to be very fast...
//
// A faster implementation that does a lot less work at runtime is definitely
// possible (it would just involve more code).


// retrieves a pointer to the field of my_chemistry that's named ``name``.
//
// This returns a NULL pointer if: my_chemistry is NULL, the field doesn't
// exist, or the field doesn't have the specified type
static void* _get_field_ptr(chemistry_data* my_chemistry, const char* name,
                            const param_list my_param_l)
{
  if (my_chemistry == NULL) { return NULL; }

  for (size_t param_index = 0; param_index < my_param_l.len; param_index++){
    const param_entry* entry = my_param_l.entries + param_index;
    if (strcmp(entry->name, name) == 0) {
      return (void*)( (char*)my_chemistry + entry->offset );
    }
  }
  return NULL;
}

int* local_chemistry_data_access_int(chemistry_data* my_chemistry,
                                     const char* param_name)
{ return (int*)_get_field_ptr(my_chemistry, param_name, _int_param_l); }

double* local_chemistry_data_access_double(chemistry_data* my_chemistry,
                                           const char* param_name)
{ return (double*)_get_field_ptr(my_chemistry, param_name, _double_param_l); }

char** local_chemistry_data_access_string(chemistry_data* my_chemistry,
                                          const char* param_name)
{ return (char**)_get_field_ptr(my_chemistry, param_name, _string_param_l); }

// define functions for accessing the names of chemistry_data:
//   param_name_int, param_name_double, param_name_string
//
// These are primarily needed for testing purposes and can be used for
// serialization.

// returns the name of the ``i``th parameter of the specified type. This returns
// NULL when there are ``i`` or fewer parameters of the specified type
static const char* _param_name(size_t i, const param_list my_param_list)
{ return (i < my_param_list.len) ? my_param_list.entries[i].name : NULL; }

const char* param_name_int(unsigned int i)
{return _param_name(i,_int_param_l);}
const char* param_name_double(unsigned int i)
{return _param_name(i,_double_param_l);}
const char* param_name_string(unsigned int i)
{return _param_name(i,_string_param_l);}

// function to query the number of parameters of chemistry_data of a given type
unsigned int grackle_num_params(const char* type_name){
  if (strcmp(type_name, "int") == 0) {
    return (unsigned int)_int_param_l.len;
  } else if (strcmp(type_name, "double") == 0) {
    return (unsigned int) _double_param_l.len;
  } else if (strcmp(type_name, "string") == 0) {
    return (unsigned int) _string_param_l.len;
  }
  return 0;
}
