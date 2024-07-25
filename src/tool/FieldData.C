#include <map>
#include "FieldData.h"
#include "grackle.h"

namespace { // anonymous namespace is local to this file

void init_grid_props(grackle_field_data* fields, int rank, int* shape) {
  fields->grid_rank = rank;
  fields->grid_dimension = new int[3];
  fields->grid_start = new int[3];
  fields->grid_end = new int[3];

  for (int i = 0; i < rank; i++) {
    fields->grid_dimension[i] = shape[i];
    fields->grid_start[i] = 0;
    fields->grid_end[i] = shape[i]-1;
  }
}

// in the future, it would probably be useful to encode the requirements for
// using different fields directly within Grackle. I have reimplemented this
// logic a handful of times...
void allocate_fields(grackle_field_data* fields, const chemistry_data& my_chem){
  std::size_t field_size = 1;
  for (int i = 0; i < fields->grid_rank; i++){
    field_size *= std::size_t(fields->grid_dimension[i]);
  }

  // skipping over a number of other choices...

  fields->density = new gr_float[field_size];
  fields->internal_energy = new gr_float[field_size];
  fields->x_velocity = nullptr;
  fields->y_velocity = nullptr;
  fields->z_velocity = nullptr;

  if (my_chem.metal_cooling == 1){
    fields->metal_density = new gr_float[field_size];
  }

  if (my_chem.primordial_chemistry >= 1) {
    fields->HI_density      = new gr_float[field_size];
    fields->HII_density     = new gr_float[field_size];
    fields->HeI_density     = new gr_float[field_size];
    fields->HeII_density    = new gr_float[field_size];
    fields->HeIII_density   = new gr_float[field_size];
    fields->e_density       = new gr_float[field_size];
  }

  if (my_chem.primordial_chemistry >= 2) {
    fields->HM_density      = new gr_float[field_size];
    fields->H2I_density     = new gr_float[field_size];
    fields->H2II_density    = new gr_float[field_size];
  }

  if (my_chem.primordial_chemistry >= 3) {
    fields->DI_density      = new gr_float[field_size];
    fields->DII_density     = new gr_float[field_size];
    fields->HDI_density     = new gr_float[field_size];
  }
}

} // anonymous namespace

grackle_field_data* impl::allocate_and_init_gr_field_data(
    const chemistry_data& my_chem, int rank, int* shape){
  grackle_field_data* fields = new grackle_field_data;
  gr_initialize_field_data(fields);
  init_grid_props(fields, rank, shape);
  allocate_fields(fields, my_chem);
  return fields;
}


// Here, we use X-Macros to create a global constant map that maps between
// fluid-field names and the pointer-to-member for fluid-fields in the
// grackle_field_data members
//
// -> a pointer-to-member is the C++ type-safe equivalent of using offsets to
//    access entries in a C struct
//
// -> DON'T IMPLEMENT THIS IN A SIMULATION CODE!!! WE ARE RELYING UPON PRIVATE
//    IMPLEMENTATION DETAILS OF GRACKLE
//
// -> If somebody wants this functionality for their simulation code, we can
//    add new functions to the grackle api (it would be implemented in C in a
//    manner more similar to the dynamic API).
//    -> If we do that, we will swap out this implementation to leverage that 
//       new API.
//    -> Adding that functionality to the grackle API would be useful for
//       supporting a stable ABI 
//
// -> Unless we're trying to support a stable ABI, this functionality is only
//    useful if you want to support arbitrary Grackle configurations in a
//    constext where the grackle_field_data struct effectively owns all the
//    field data.
//    -> The usage of field data in most simulation codes is quite a bit more
//       constrained. So this probably isn't very useful

#define MK_PAIR(FIELD) {#FIELD, &grackle_field_data::FIELD},

static const std::map<std::string,
                      gr_float* grackle_field_data::*> field_mapping_ {
  #define ENTRY(MEMBER_NAME) MK_PAIR(MEMBER_NAME)
  #include "grackle_field_data_fdatamembers.def"
  #undef ENTRY
};

void destroy_selfcontained_field_data(grackle_field_data* my_fields) {
  delete[] my_fields->grid_dimension;
  delete[] my_fields->grid_start;
  delete[] my_fields->grid_end;

  for (auto kv_pair: field_mapping_){
    gr_float* grackle_field_data::* ptr_to_mem = kv_pair.second;
    gr_float* ptr = my_fields->*ptr_to_mem;
    if (ptr != nullptr) { delete[] ptr; }
  }

  delete my_fields;
}


void clone_field_data(grackle_field_data* dest, grackle_field_data* src) {
  GRCLI_REQUIRE((dest != nullptr) && (src != nullptr), "args can't be null");
  GRCLI_REQUIRE(dest->grid_rank == src->grid_rank, "the ranks must match");

  std::size_t size_src = 1;
  std::size_t size_dest = 1;
  for( int i = 0; i < dest->grid_rank; i++) {
    size_src *= src->grid_dimension[i];
    size_dest *= dest->grid_dimension[i];
  }
  GRCLI_REQUIRE(size_src == size_dest,
                "the number of allocated elements doesn't match");

  for( int i = 0; i < dest->grid_rank; i++) {
    dest->grid_dimension[i] = src->grid_dimension[i];
    dest->grid_start[i] = src->grid_start[i];
    dest->grid_end[i] = src->grid_end[i];
  }


  for (auto kv_pair: field_mapping_){
    gr_float* grackle_field_data::* ptr_to_mem = kv_pair.second;
    const gr_float* src_ptr = src->*ptr_to_mem;
    gr_float* dest_ptr = dest->*ptr_to_mem;
    if ((src_ptr == nullptr) && (dest_ptr == nullptr)) {
      continue;
    } else if ((src_ptr != nullptr) && (dest_ptr != nullptr)) {
      std::memcpy(dest_ptr, src_ptr, size_src);
    } else {
      GRCLI_ERROR("The src or dest object has NULL pointer for the %s while "
                  "the other has a non-NULL pointer", kv_pair.first.c_str());
    }
  }
}

View<gr_float> FieldData::view(std::string_view field_name) const {
    GRCLI_ERROR("NOT IMPLEMENTED YET!!!!");
}



