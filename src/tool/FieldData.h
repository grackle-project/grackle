#ifndef FIELD_DATA_H
#define FIELD_DATA_H

#include <array>
#include <string_view>
#include <vector>


#include "view.h"

#include "grackle.h"

namespace impl {

/// allocate a grackle_field_data struct tha owns its own allocations
grackle_field_data* allocate_and_init_gr_field_data(
    const chemistry_data& my_chem, int rank, int* shape);

/// deallocate all memory associated with a grackle_field_data struct that owns
/// all associated memory
void destroy_selfcontained_field_data(grackle_field_data* my_fields);

} // impl namespace

/// @brief wraps grackle_field_data and all underlying allocations
///
/// In many simulation codes, a ``grackle_field_data`` instance doesn't own the
/// underlying allocations of the field data. But for this application, that is
/// necessary
class FieldData {
  grackle_field_data* grackle_fields_;

public:

  // remove the default constructor
  FieldData() = delete;

  // remove the copy-constructor and copy-assignment
  FieldData(const FieldData&) = delete;
  FieldData& operator=(const FieldData&) = delete;

  // for concreteness, be clear that move-constructor and move-assignment work
  FieldData(FieldData&&) = default;
  FieldData& operator=(FieldData&&) = default;

  FieldData(std::array<int, 3> shape, const chemistry_data& my_chem)
    : grackle_fields_{impl::allocate_and_init_gr_field_data(my_chem, 3,
                                                            shape.data())}
  {}

  ~FieldData()
  { impl::destroy_selfcontained_field_data(grackle_fields_); }

  View<gr_float> view(std::string_view field_name) const;

  std::array<int, 3> grid_dimensions() const {
    return {grackle_fields_->grid_dimension[0],
            grackle_fields_->grid_dimension[1],
            grackle_fields_->grid_dimension[2]};
  }

  grackle_field_data* get_ptr() noexcept {return grackle_fields_; }
  const grackle_field_data* get_ptr() const noexcept { return grackle_fields_; }

};

#endif /* FIELD_DATA_H */
