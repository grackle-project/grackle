#ifndef CHEMISTRY_DATA_H
#define CHEMISTRY_DATA_H

#include <memory> // std::unique_ptr
#include <string>
#include <string_view>
#include <type_traits> // std::is_same
#include <unordered_set>
#include <utility> // std::pair
#include <vector>

#include "grackle.h"

#include "utils.h"

/// Wrapper around the chemistry_data struct
///
/// @note
/// This is a stripped-down (and slightly updated) version of Enzo-E's
/// `GrackleChemistryData` class
///
/// @note
/// Instances of this class satisfy 2 invariants:
///   1. the stored pointer to the ``chemistry_data`` struct is always valid;
///      it's not allowed to be a ``nullptr``. While lazily initializing the
///      pointer would make the move constructor cheaper, it's probably not a
///      worthwhile tradeoff (it would be REALLY easy to forget to lazily
///      initialize the pointer and then cause a bug)
///   2. the str_allocs_ attribute manages the lifetime of all strings
///      stored in fields of the ``chemistry_data`` struct. The only
///      exception is the default value (which may be a string literal or a
///      ``nullptr``)
class ChemistryData {

  /// pointer to the chemistry_data struct
  std::unique_ptr<chemistry_data> ptr_;

  /// manages the lifetime of string fields used in chemistry_data
  ///
  /// We explicitly avoid using std::string because SSO (small string
  /// optimization) could cause all sorts of bugs in the future that would be
  /// VERY hard to debug
  ///
  /// The number of entries in this parameter is currently small. In the
  /// future, it may be sensible to replace vector with unordered_set if the 
  /// number grows
  std::vector<std::unique_ptr<char[]>> str_allocs_;

  bool try_set_int_(const std::string& field, int value) {
    int* field_ptr = local_chemistry_data_access_int(ptr_.get(), field.c_str());
    if (!field_ptr) { return false; }
    *field_ptr = value;
    return true;
  }


  bool try_set_dbl_(const std::string& field, double value) {
    double* field_ptr = local_chemistry_data_access_double(ptr_.get(),
                                                           field.c_str());
    if (!field_ptr) { return false; }
    *field_ptr = value;
    return true;
  }

  inline bool try_set_str_(const std::string& field, std::string_view value);

public:
  /// construct a new ChemistryData
  ChemistryData()
    : ptr_{new chemistry_data}, str_allocs_()
  { local_initialize_chemistry_parameters(ptr_.get()); }

  /// destroy a ChemistryData instance
  ~ChemistryData() = default;

  // eliminate the ability to perform deepcopies (for now)
  ChemistryData(const ChemistryData& other) = delete;
  ChemistryData& operator= (const ChemistryData&) = delete;

  /// move constructor
  ChemistryData(ChemistryData&& other)
    : ChemistryData()
  { this->swap(other); }

  /// move-assignment operation
  ChemistryData& operator= (ChemistryData&& other) noexcept
  { this->swap(other); return *this; }

  /// exchange the contents of this with other
  void swap(ChemistryData& other) noexcept
  { ptr_.swap(other.ptr_); str_allocs_.swap(other.str_allocs_); }

  /// returns a pointer to the managed chemistry_data struct
  ///
  /// This is primarily intended to be used when calling a grackle function.
  /// The pointer is only valid during the ChemistryData's lifetime
  inline chemistry_data* get_ptr() { return ptr_.get(); }
  inline const chemistry_data* get_ptr() const { return ptr_.get(); }

  /// try to update a parameter value stored in the chemistry_data struct.
  ///
  /// @retval true the stored value was succesfully updated
  /// @retval false there is no know parameter of the specified type
  ///
  /// @note
  /// specializations of this method follow the class declaration
  template<class T>
  bool try_set(const std::string& field, T value) {
    using cleanT = std::remove_cv_t<std::remove_reference_t<T>>;

    if constexpr (std::is_same_v<cleanT,int>) {
      return try_set_int_(field, value);
    } else if constexpr (std::is_same_v<cleanT,double>) {
      return try_set_dbl_(field, value);
    } else if constexpr (std::is_same_v<cleanT,std::string>) {
      return try_set_str_(field, value);
    } else if constexpr (std::is_same_v<cleanT,std::string_view>) {
      return try_set_str_(field, value);
    } else if constexpr (std::is_same_v<cleanT,char*>) {
      return try_set_str_(field, value);
    } else {
      GRCLI_ERROR("the value must be an int, a double, or a string");
    }
  }

};


inline bool ChemistryData::try_set_str_(const std::string& field,
                                        std::string_view value)
{

  // NOTE: we should NOT directly modify characters held by field_ptr
  char ** field_ptr = local_chemistry_data_access_string(ptr_.get(),
                                                         field.c_str());

  if (!field_ptr) { return false; }

  // deallocate the existing value (if applicable)
  if ((*field_ptr) != nullptr){

    // check whether *field_ptr matches the address of any pointers held within
    // str_allocs_, if so delete that pointer from str_allocs_
    for (std::size_t i = 0; i < str_allocs_.size(); i++){

      if (str_allocs_[i].get() == *field_ptr){
        // we will only enter this part of the loop if *field_ptr doesn't refer
        // to a string literal (that could have been set as a default value)
        str_allocs_.erase(str_allocs_.begin() + i);
        break;
      }
    }

  }

  // allocate a new c-string and copy data from value into it
  const std::size_t length_with_nul = value.size() + 1;
  std::unique_ptr<char[]> new_alloc(new char[length_with_nul]);
  std::memcpy(new_alloc.get(), value.data(), length_with_nul-1);
  new_alloc.get()[length_with_nul-1] = '\0';

  // update the field of the chemistry_data struct
  (*field_ptr) = new_alloc.get();

  // finally, move the newly allocated c-string into str_allocs_
  str_allocs_.push_back(std::move(new_alloc));
  return true;
}


#endif /* CHEMISTRY_DATA_H */
