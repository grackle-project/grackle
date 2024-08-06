#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <regex>
#include <string>
#include <string_view>

#include "handle_gr_parameters.h"
#include "utils.h"

namespace { // stuff inside anonymous namespace is local to this file

// the white-space on either side of the equal-sign is ONLY meaningful if the
// user properly quotes the key-value pair
const std::regex kv_regex("^([a-zA-Z_0-9]+)[\\s]*=[\\s]*(.*)$");

const std::regex flt_regex("^" FLT_PATTERN "$");
const std::regex int_regex("^\\d+$");
const std::regex str_regex("^" STRING_PATTERN "$");


// the returned string is empty if there aren't any problems
template<typename Fn>
std::string for_each_char_(std::string_view s, Fn fn) {
  // we assume that s[0] == '"' and s[s.size()-1] == '"'
  int i = 1;
  int stop = s.size()-1;
  while (i < stop) {
    if (s[i] == '\\'){
      if ((i+1) == stop) return {"lone backslash character"};

      switch (s[i+1]) {
        case '"':  { fn('"'); break; }
        case '\\': { fn('\\'); break; }
        case 'n':  { fn('\n'); break; }
        case 't':  { fn('\t'); break; }
        default:
          return (std::string("invalid escape-sequence: \\") + s[i+1]);
      }
      i+=2;
    } else {
      fn(s[i]);
      i++;
    }
  }
  return {};
}

/// returns true for success
///
/// when this returns false and rslt is empty, then s simply didn't correspond
/// to a string. When this returns false and rslt isn't empty, then rslt
/// specifies an error.
bool parse_string_val_(std::string s, std::string& rslt){
  if ((s.size() > 2) && !std::regex_match(s, str_regex)) {
    rslt.clear();
    return false;
  }

  // count the number of characters in the string
  int count = 0;
  std::string err_msg = for_each_char_(s, [&count](char){ count++; });

  if (!err_msg.empty()) {
    rslt = err_msg;
    return false;
  } else {
    rslt.clear();
    rslt.reserve(count);
    for_each_char_(s, [&rslt](char ch){ rslt.push_back(ch); });
    return true;
  }
}

} // anonymous namespace


CliParamSpec::KVPair CliParamSpec::parse_param_(std::string_view token){

  // copy the token
  std::string token_str(token);
  //std::printf("parsing %s\n", token_str.c_str());

  CliParamSpec::KVPair out;
  std::smatch kv_match;
  if (!std::regex_match(token_str, kv_match, kv_regex)) {
    std::fprintf(stderr,
        "grackle-parameters should have the form `<key>=<value>`. The "
        "`%s` argument does not meet these expectations.", token_str.c_str());
    std::exit(1);
  }

  // match 0 is full string
  // extract the key-name
  std::string_view key(token.data() + kv_match.position(1), kv_match.length(1));
  
  // extract the substring corresponding to the value and then parse it
  std::string value_str = kv_match[2].str();

  value_variant value;
  if (std::regex_match(value_str, int_regex)) {
    int tmp_val = std::stoi(value_str);
    value = tmp_val;
  } else if (std::regex_match(value_str, flt_regex)) {
    double tmp_val = std::stod(value_str);
    value = tmp_val;
  } else {
    std::string tmp_val;
    bool is_str = parse_string_val_(value_str, tmp_val);

    if (is_str) {
      value = tmp_val;
    } else if (!tmp_val.empty()) {
      // in this case, tmp_val encodes an error message
      std::fprintf(stderr,
          "An error occured while parsing the right-hand-side of the `%s` "
          "grackle-parameter argument as a string: %s\n",
          token_str.c_str(), tmp_val.c_str());
      std::exit(1);
    } else {
      std::string key_str{key};
      std::fprintf(stderr,
          "The expression on the right-hand-side of the `%s`\n"
          "grackle-parameter argument is invalid.\n"
          "-> It must be an integer, floating-point or doubly-quoted str\n"
          "-> If you are using bash and the argument looked like\n"
          "       %s=\"%s\"\n"
          "   then there's a problem with quote-escaping. You may want to "
          "try passing:\n"
          "       '%s=\"%s\"'\n",
          token_str.c_str(), key_str.c_str(), value_str.c_str(),
          key_str.c_str(), value_str.c_str());
      std::exit(1);
    }
  }

  return {key, value};
}

static const char* fetch_gr_parameter_type(std::string_view key) {
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
static void init_gr_params_(const CliParamSpec& param_spec,
                            ChemistryData& my_chem)
{

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

FullGrackleSolverPack create_full_grackle_solver(
  std::optional<CliParamSpec> maybe_parameter_spec)
{
  // it may be beneficial to introduce the idea of presets to avoid needing to
  // initialize parameters like use_grackle

  ChemistryData my_chem;
  if (maybe_parameter_spec.has_value()) {
    init_gr_params_(maybe_parameter_spec.value(), my_chem);
  }

  // in the future, we will support customization of units!
  double initial_redshift = 1.0;
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 3.0857e18; // ~ 1 parsec
  my_units.time_units = 3.15e13; // ~ 1 Myr
  my_units.a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;

  return FullGrackleSolverPack(std::move(my_chem), my_units);
}
