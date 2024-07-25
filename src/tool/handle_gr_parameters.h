#ifndef TOOL_GR_PARAMETERS_H
#define TOOL_GR_PARAMETERS_H

#include <string_view>
#include <string>
#include <variant>

#include "grackle.h"

#include "ChemistryData.h"
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
    for(const char* const* ptr = this->begin + 1; ptr != end; ++ptr) {
      KVPair kv_pair = this->parse_param_(*ptr);
      f(kv_pair.key, kv_pair.val);
    }
  }

};

#define LEADING_PARAM_ARG "--par-start"

inline bool check_start_cli_paramspec(const char* arg) {
  return std::string_view(LEADING_PARAM_ARG) == arg;
}

inline CliParamSpec parse_cli_paramspec(const char* const* begin,
                                        const char* const* end){
  GRCLI_REQUIRE(check_start_cli_paramspec(*begin),
                "The first argument must be \"%s\"", LEADING_PARAM_ARG);

  const char* const* tmp = begin;
  while(true){
    tmp++;

    if ((tmp == end) || (std::string_view("--par-stop") == *tmp)) {
      if (tmp == (begin + 1)){
        fprintf(stderr, "the \"%s\" flag doesn't start a group of parameters\n",
                LEADING_PARAM_ARG);
        std::exit(1);
      }
      return {begin, tmp};
    }
  }

}


inline const char* fetch_gr_parameter_type(std::string_view key) {
  using F = const char*(unsigned int);
  F* fn_list[3] = {param_name_int, param_name_double, param_name_string};
  const char* type_list[3] = {"int", "double", "string"};

  for (int i = 0; i < 3; i++) {
    F* fn = fn_list[i];
    const char* type_name = type_list[i];
    unsigned int num_pars = grackle_num_params(type_name);

    for (unsigned int j = 0; j < num_pars; j++) {
      if (fn(j) == key) return type_name;
    }
  }
  return nullptr;
}

// the plan would be to eventually support a hybrid approach (e.g. we could use
// CliParamSpec and dumped parameters
inline void init_gr_params(const CliParamSpec& param_spec,
                           ChemistryData& my_chem)
{
  /*
  auto fn = [](std::string_view key, value_variant value) {
    std::string key_str(key.data(), key.size());
    if (std::holds_alternative<int>(value)) {
      std::printf("-> %s: %d (INTEGER)\n", key_str.c_str(), std::get<int>(value));
    } else if (std::holds_alternative<double>(value)) {
      std::printf("-> %s: %g (DOUBLE)\n", key_str.c_str(), std::get<double>(value));
    } else if (std::holds_alternative<std::string>(value)) {
      std::printf("-> %s: '%s' (STRING)\n", key_str.c_str(),
                  std::get<std::string>(value).c_str());
    } else {
      GRCLI_ERROR("SOMETHING IS VERY_WRONG!!!!");
    }
  };
  */


  auto fn = [&my_chem](std::string_view key, value_variant value) {
    std::string key_str(key.data(), key.size());
    bool success;
    const char* specified_type = nullptr;
    if (std::holds_alternative<int>(value)) {
      success = my_chem.try_set(key_str, std::get<int>(value));
      specified_type = "int";
    } else if (std::holds_alternative<double>(value)) {
      success = my_chem.try_set(key_str, std::get<double>(value));
      specified_type = "double";
    } else if (std::holds_alternative<std::string>(value)) {
      success = my_chem.try_set(key_str, std::get<std::string>(value));
      specified_type = "string";
    } else {
      GRCLI_ERROR("Something is very wrong! Encountered unknown type");
    }


    if (!success) {
      const char* expected_type = fetch_gr_parameter_type(key);
      if (expected_type == nullptr) {
        std::fprintf(stderr, "\"%s\" is not a known grackle parameter\n",
                    key_str.c_str());
      } else {
        std::fprintf(stderr, 
            "the \"%s\" grackle parameter was specified with \"%s\" value. "
            "The value should be of type \"%s\"\n",
            key_str.c_str(), specified_type, expected_type);
      }
      std::exit(1);
    }
  };

  param_spec.for_each(fn);

}


#endif /* TOOL_GR_PARAMETERS_H */
