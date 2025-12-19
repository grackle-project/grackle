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

#include <cstddef>
#include <cstring>
#include "grackle_chemistry_data.h"
#include "status_reporting.h"

namespace {  // stuff in an anonymous namespace is only visible to this file

/// This is collection of enumerators, with an enumerator named for EVERY
/// parameter. Importantly, each enumerator has a unique value (that we use
/// in a giant switch-statement).
enum ParamEnum {
#define ENTRY(FIELD, TYPE, DEFAULT_VAL) ParamEnum_ ## FIELD,
#include "grackle_chemistry_data_fields.def"
#undef ENTRY

  ParamEnum_NUM_ENTRIES  // <- always last (so it specifies number of entries)
};

// initialize int_param_l_, double_param_l_ & string_param_l_. They are sized
// lists that respectively hold entries for each int, double, and string field
// in chemistry_data. These are only used internally
struct param_entry {
  const enum ParamEnum enum_val;
  const char * name;
};

struct param_list {
  const std::size_t len;
  const param_entry * entries;
};

#define INIT_PAR_LIST(ENTRIES) {sizeof(ENTRIES) / sizeof(param_entry), ENTRIES}

#define LIST_INT_INT(FIELD) {ParamEnum_ ## FIELD, #FIELD},
#define LIST_INT_DOUBLE(FIELD) /* ... */
#define LIST_INT_STRING(FIELD) /* ... */

#define LIST_DOUBLE_INT(FIELD) /* ... */
#define LIST_DOUBLE_DOUBLE(FIELD) {ParamEnum_ ## FIELD, #FIELD},
#define LIST_DOUBLE_STRING(FIELD) /* ... */

#define LIST_STRING_INT(FIELD) /* ... */
#define LIST_STRING_DOUBLE(FIELD) /* ... */
#define LIST_STRING_STRING(FIELD) {ParamEnum_ ## FIELD, #FIELD},

constexpr param_entry int_param_entries_[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) LIST_INT_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
constexpr param_list int_param_l_ = INIT_PAR_LIST(int_param_entries_);

constexpr param_entry double_param_entries_[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) LIST_DOUBLE_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
constexpr param_list double_param_l_ =
  INIT_PAR_LIST(double_param_entries_);

constexpr param_entry string_param_entries_[] = {
  #define ENTRY(FIELD, TYPE, DEFAULT_VAL) LIST_STRING_ ## TYPE(FIELD)
  #include "grackle_chemistry_data_fields.def"
  #undef ENTRY
};
constexpr param_list string_param_l_ = INIT_PAR_LIST(string_param_entries_);


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


/// retrieves a pointer to the field of my_chemistry that's named ``name``.
///
/// This returns a nullptr pointer if: my_chemistry is nullptr, the field doesn't
/// exist, or the field doesn't have the specified type
static void* get_field_ptr_(chemistry_data* my_chemistry, const char* name,
                            const param_list my_param_l)
{
  if (my_chemistry == nullptr) { return nullptr; }

  // lookup the enum value associated with name
  enum ParamEnum enum_val = ParamEnum_NUM_ENTRIES;
  for (std::size_t param_index = 0; param_index < my_param_l.len; param_index++){
    const param_entry* entry = my_param_l.entries + param_index;
    if (std::strcmp(entry->name, name) == 0) {
      enum_val = entry->enum_val;
      break;
    }
  }

  // now use a switch statement to actually return the appropriate pointer
  switch (enum_val) {
    #define ENTRY(FIELD, TYPE, DEFAULT_VAL)                              \
      case ParamEnum_ ## FIELD: { return (void*)(&my_chemistry->FIELD); }
    #include "grackle_chemistry_data_fields.def"
    #undef ENTRY

    case ParamEnum_NUM_ENTRIES: return nullptr;  // <- unknown name
  }
  GR_INTERNAL_UNREACHABLE_ERROR();
}

}  // anonymous namespace

extern "C" {

int* local_chemistry_data_access_int(chemistry_data* my_chemistry,
                                     const char* param_name)
{ return (int*)get_field_ptr_(my_chemistry, param_name, int_param_l_); }

double* local_chemistry_data_access_double(chemistry_data* my_chemistry,
                                           const char* param_name)
{ return (double*)get_field_ptr_(my_chemistry, param_name, double_param_l_); }

char** local_chemistry_data_access_string(chemistry_data* my_chemistry,
                                          const char* param_name)
{ return (char**)get_field_ptr_(my_chemistry, param_name, string_param_l_); }

} // extern "C"

// define functions for accessing the names of chemistry_data:
//   param_name_int, param_name_double, param_name_string
//
// These are primarily needed for testing purposes and can be used for
// serialization.

// returns the name of the ``i``th parameter of the specified type. This returns
// nullptr when there are ``i`` or fewer parameters of the specified type
static const char* param_name_(std::size_t i, const param_list my_param_list)
{ return (i < my_param_list.len) ? my_param_list.entries[i].name : nullptr; }

extern "C" {

const char* param_name_int(unsigned int i)
{return param_name_(i,int_param_l_);}
const char* param_name_double(unsigned int i)
{return param_name_(i,double_param_l_);}
const char* param_name_string(unsigned int i)
{return param_name_(i,string_param_l_);}

// function to query the number of parameters of chemistry_data of a given type
unsigned int grackle_num_params(const char* type_name){
  if (type_name == nullptr) {
    return 0;
  } else if (std::strcmp(type_name, "int") == 0) {
    return (unsigned int) int_param_l_.len;
  } else if (std::strcmp(type_name, "double") == 0) {
    return (unsigned int) double_param_l_.len;
  } else if (std::strcmp(type_name, "string") == 0) {
    return (unsigned int) string_param_l_.len;
  }
  return 0;
}

} // extern "C"
