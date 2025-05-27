/// @file misc.hpp
/// @brief implement miscellaneous types/functions

#ifndef GRTEST_EVOLVE_MISC_HPP
#define GRTEST_EVOLVE_MISC_HPP

#include <grackle.h>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>  // std::pair
#include <vector>

#include "evolve.hpp"

namespace grtest::impl {

int elements_per_field_ptr(const grackle_field_data& my_fields) {
  int size = 1;
  for (int i = 0; i < my_fields.grid_rank; i++) {
    size *= my_fields.grid_dimension[i];
  }
  return size;
}

bool has_ghost_padding(const grackle_field_data& my_fields) {
  for (int i = 0; i < my_fields.grid_rank; i++) {
    if ((my_fields.grid_start[i] != 0) ||
        (my_fields.grid_end[i] + 1 != my_fields.grid_dimension[i])) {
      return true;
    }
  }
  return false;
}

/// creates a function based on the minimum cooling time
///
/// @param safety_factor
/// @param pack Holds the grackle types
std::function<double()> make_tcool_dt_func(
  double safety_factor, grtest::GrackleTypePack pack
) {
  const gr_float max_tcool = std::numeric_limits<gr_float>::max();

  // initialize a vector to hold intermediate values
  int size = grtest::impl::elements_per_field_ptr(*pack.my_fields);
  std::vector<gr_float> vec(static_cast<std::size_t>(size), max_tcool);

  // define an inline function
  auto fn = [vec = std::move(vec), safety_factor, pack, max_tcool]() mutable {
    // calculate the cooling time
    int rv = local_calculate_cooling_time(
      pack.my_chem, pack.my_rates, pack.my_units, pack.my_fields, vec.data()
    );
    if (rv != GR_SUCCESS) { return -1.0; }

    // find the minimum value
    gr_float minimum = max_tcool;
    for (gr_float val : vec) {
      gr_float abs_val = std::fabs(val);
      minimum = (abs_val < minimum) ? abs_val : minimum;
    }

    return safety_factor * minimum;
  };

  return std::function<double()>(fn);
}

/// A callable that indicates whether we have triggered the stopping condition
class FieldStopTrigger {
  FieldStopThresh kind;
  double threshold_cgs;
  int n_elements;
  std::vector<gr_float> buffer; // used for temperature
  grtest::GrackleTypePack pack;

  FieldStopTrigger() = default;
public:

  /// Returns true if we triggered the stopping condition, false if we haven't
  /// and an empty optional to denote an error
  std::optional<bool> operator()() {
    if (kind == FieldStopThresh::None) {
      return {false};
    } else if (kind == FieldStopThresh::MinTemperature) {
      int rv = local_calculate_temperature(
        pack.my_chem, pack.my_rates, pack.my_units, pack.my_fields,
        buffer.data()
      );
      if (rv != GR_SUCCESS) { return {}; }
      bool out = false;
      for (int i = 0; i < n_elements; i++) {
        if (buffer[i] <= threshold_cgs) { out = true; }
      }
      return {out};
    } else if (kind == FieldStopThresh::MinTemperature) {
      double density_units = pack.my_units->density_units;
      gr_float* ptr = pack.my_fields->density;
      bool out = false;
      for (int i = 0; i < n_elements; i++) {
        if (ptr[i] * density_units >= threshold_cgs) { out = true; }
      }
      return {out};
    } else {
      return {};
    }
  }

  /// create a new FieldStopTrigger
  static std::pair<std::string, FieldStopTrigger> create(
    FieldStopThresh kind, double threshold_cgs, grtest::GrackleTypePack pack
  ) {
    FieldStopTrigger trigger;
    if (kind != FieldStopThresh::None) {
      if (has_ghost_padding(*pack.my_fields)) {
        return {"FieldStopTrigger doesn't support ghost padding yet", trigger};
      } else if (threshold_cgs <= 0) {
        return {"FieldStopTrigger needs positive threshold", trigger};
      }
    }
    trigger.kind = kind;
    trigger.threshold_cgs = threshold_cgs;
    trigger.n_elements = elements_per_field_ptr(*pack.my_fields);
    trigger.buffer = std::vector<gr_float>((std::size_t)trigger.n_elements);
    trigger.pack = pack;
    return {"", trigger};
  }
};

}  // namespace grtest::impl

#endif /* GRTEST_EVOLVE_STOPTRIGGER_HPP */

