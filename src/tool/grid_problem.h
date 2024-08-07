#ifndef GRID_PROBLEM_H
#define GRID_PROBLEM_H

#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "grackle.h"

#include "utils.h"



namespace valseq{

struct Single{
  double lo, hi;
  int num;
  Single (double val) : lo(val), hi(val), num(1) {}
  double get(int i) const { return lo; }
};

// todo: refactor to support repeating an array of values
struct Repeat {
  double lo, hi;
  int num;
  Repeat (double val, int num) : lo(val), hi(val), num(num)
  {
    GRCLI_REQUIRE(num > 0, "num must exceed 0");
  }
  double get(int i) const { return lo; }
};

struct Geomspace {
  double lo, hi;
  int num;

private:
  double loglo, loghi;
  double log_delta;

public:
  Geomspace(double lo, double hi, int num)
    : lo(lo), hi(hi), num(num)
  {
    GRCLI_REQUIRE(lo > 0, "lo must be positive");
    GRCLI_REQUIRE(lo < hi, "hi, %g, doesn't exceed lo, %g", hi, lo);
    GRCLI_REQUIRE(num > 1, "num must exceed 1");

    loglo = std::log(lo);
    loghi = std::log(hi);
    log_delta = (loghi - loglo) / (num - 1);
  }

  double get(int i) const {
    if (i == 0) {
      return lo;
    } else if (i+1 == num) {
      return hi;
    } else {
      return std::exp(loglo + i * log_delta);
    }
  }

};

}  // namespace valseq

using seq_variant = std::variant<valseq::Single, valseq::Geomspace>;

// A scenario is basically the recipe for initializing the field data
// -> at the moment, the grid scenario is the only kind. It uses a sequence
//    of target values along each axis to initialize the grid
// -> in principle, we might support other kinds of scenarios

namespace scenario{

enum class QuantityKind { temperature, mass_density, metallicity };

struct AxisSpec {
  QuantityKind quantity;
  seq_variant seq;
};

// not currently used (but could be specified in the future)
// enum class IonizationAssumption {ionized, neutral, fixed_mmw };

struct CliGridSpec {
  /// records whether any component is specifed
  bool any_specified;

  /// records the quantity sequence along a specified axis
  std::optional<AxisSpec> ax_props[3];
};


/// Tries to parse arg to see if it corresponds to a component of
/// GridScenarioSpec
///
/// * If so, this will (i) parse the argument, (ii) modify the `grid_spec`
///   function-arg accordingly, and (iii) return true
///
/// * Otherwise, this returns false
///
/// * If the argument is invalid, the program will abort with an error
bool try_parse_cli_grid_component(const char* arg,
                                  CliGridSpec& grid_spec);


FieldData initialize_grid(const chemistry_data& my_chem,
                          const code_units& initial_units,
                          const CliGridSpec& scenario);

} // namespace scenario

//std::optional<int> required_primordial_chemistry() const {
//  return {};
//}

#endif /* GRID_PROBLEM_H */
