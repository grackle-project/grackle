#include <optional>
#include <regex>
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


/// Parse the a sequence command-line argument
inline seq_variant get_seq(std::string s){

  const std::regex flt_regex("^" FLT_PATTERN "$");

  if (std::regex_match(s, flt_regex)) {
    double val = std::stod(s);
    valseq::Single seq(val);
    return { seq };

  } else if (starts_with(s, "geomspace")) { // "geomspace;lo;hi;num"
    const char* pattern = (
        "^([a-zA-Z]+);(" FLT_PATTERN ");(" FLT_PATTERN ");([-+]?\\d+)$");
    std::regex geomspace_regex(pattern);
    std::smatch match;
    if (std::regex_match(s, match, geomspace_regex)) {
      // match 0 is full string
      // match 1 is "geomspace"

      // match 2 is the lo string
      double lo = std::stod(match[2].str());
      // matches 3-5 are ALWAYS various parts of float (even if we the string
      // doesn't have characters that match every subgroup)

      // match 6 is the hi string
      double hi = std::stod(match[6].str());
      // matches 7-9 are ALWAYS various parts of a float

      // match 10 is always the number of entries
      int num = std::stoi(match[10].str());
      printf("geomspace;%g;%g;%d\n", lo,hi,num);
      valseq::Geomspace seq(lo,hi,num);
      return {seq};
    }

    GRCLI_ERROR("geomspace must be formatted as geomspace;{lo};{hi};{num}");
  } else {
    GRCLI_ERROR("unrecognized value");
  }
}

enum class QuantityKind { temperature, mass_density, metallicity };

struct AxisSpec {
  QuantityKind quantity;
  seq_variant seq;
};

enum class IonizationAssumption {ionized, neutral, fixed_mmw };

struct GridScenarioSpec {
  std::optional<AxisSpec> ax_props[3];
  // IonizationAssumption could be something in the future
};

grackle_field_data* initialize_grid(chemistry_data* my_chemistry,
                                    code_units* initial_units,
                                    GridScenarioSpec scenario){
  return nullptr;
}



//std::optional<int> required_primordial_chemistry() const {
//  return {};
//}

