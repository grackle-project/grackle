#ifndef TOOL_GR_PARAMETERS_H
#define TOOL_GR_PARAMETERS_H

#include <memory>
#include <optional>
#include <string_view>
#include <string>
#include <utility>
#include <variant>

#include "grackle.h"

#include "ChemistryData.h"
#include "CliParser.h"
#include "utils.h"


using value_variant = std::variant<int, double, std::string>;

/// tracks the sequence of cli arguments that specify grackle parameters
///
/// @note
/// This is intended to track pointers inside the argv argument that main
/// receives. If the order within argv, this becomes invalidated. Likewise
/// this would also become invalidated, if argv were somehow deallocated...
struct CliParamSpec{
  const char* const * begin;
  const char* const * end;

private:
  /// internal type used to hold the parsed key-value pair
  struct KVPair{ std::string_view key; value_variant val; };

  /// parse the key-value pair
  ///
  /// @note
  /// This assumes that each key-value pair was specified as a single token.
  /// This also assumes that the token's lifetime will be at least as long as
  /// the string_view in the returned object
  ///
  /// @note
  /// due to some limitations in the C++ standard library's regex
  /// implementation, it needs to operate on std::string rather than
  /// std::string_view. Internally, we also need to allocate std::string. As a
  /// result, things are quite a bit slower than they could be. 
  /// - If we decide that performance is an issue, then there are definitely a
  ///   handful of optimizations that could be made
  /// - (but it may be better to adopt an external library that can operate on
  ///   fragments).
  static KVPair parse_param_(std::string_view token);

public:
  /// calls the function on each key-value pair. The function should expect the
  /// arguments (std::string_view key, value_variant value)
  template<class PairFn>
  void for_each(PairFn f) const {
    for(const char* const* ptr = this->begin; ptr != end; ++ptr) {
      KVPair kv_pair = this->parse_param_(*ptr);
      f(kv_pair.key, kv_pair.val);
    }
  }

};


/// Checks whether the specified arg(s) correspond to a CliParamSpec:
///
/// * If so, this will parse the arguments and update the param_spec
///   argument. This function will return true
///
/// * Otherwise, this returns false
inline bool try_parse_cli_paramspec(const char* leading_arg,
                                    CliParser& parser,
                                    std::optional<CliParamSpec>& param_spec) {
  if (std::string_view(leading_arg) != "--par-start") {
    return false;
  } else if (param_spec.has_value()) {
    std::fprintf(stderr, "%s isn't allowed to be provided more than once\n",
                 leading_arg);
    std::exit(1);
  }

  const std::string_view sentinel("--par-stop");
  if ((not parser.has_next()) || (parser.peek() == sentinel)) {
    fprintf(stderr, "the \"%s\" flag doesn't start a group of parameters\n",
            leading_arg);
    std::exit(1);

  } else {
    char * const * begin = parser.next_argv_ptr();
    char * const * latest = begin;
    while (sentinel != *latest) {
      if (!parser.has_next()) {
        latest++;
        break;
      }
      latest = parser.next_argv_ptr();
    }

    CliParamSpec tmp{begin, latest};
    param_spec = std::make_optional<CliParamSpec>(std::move(tmp));
    return true;
  }
}

/// Class that holds the full configuration of the Grackle-Solver
class FullGrackleSolverPack {
  ChemistryData parameters_;
  std::unique_ptr<chemistry_data_storage> rates_;
  code_units initial_units_;

public:
  FullGrackleSolverPack() = delete;

  FullGrackleSolverPack(ChemistryData&& parameters, code_units initial_units)
    : parameters_(std::move(parameters)),
      rates_(new chemistry_data_storage),
      initial_units_(initial_units)
  {
    if (local_initialize_chemistry_data(parameters_.get_ptr(), rates_.get(),
                                        &initial_units_) == GR_FAIL) {
      GRCLI_ERROR("Error in local_initialize_chemistry_data");
    }
  }

  ~FullGrackleSolverPack() {
    local_free_chemistry_data(parameters_.get_ptr(), rates_.get());
  }

  chemistry_data* chemistry_data() { return parameters_.get_ptr(); }

  chemistry_data_storage* chemistry_storage() { return rates_.get(); }

  code_units initial_units() { return initial_units_; }
};

/// Fully initialize all parts of the grackle solver
///
/// @note
/// Realistically, this may need to be refactored if we tack-on more features
FullGrackleSolverPack create_full_grackle_solver(
  std::optional<CliParamSpec> maybe_parameter_spec);


#endif /* TOOL_GR_PARAMETERS_H */
