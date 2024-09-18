#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <regex>
#include <string_view>
#include <string>

#include "operation.h"
#include "utils.h"

namespace{ // anonymous namespace contents only visible in this file

struct OperationProps{
  OperationKind kind;
  const char* name;
  bool requires_dt;
};

const OperationProps op_prop_list[] = {
  {OperationKind::calc_cooling_time, "calc-cooling-time", false},
  {OperationKind::calc_dust_temperature, "calc-dust-temperature", false},
  {OperationKind::calc_pressure, "calc-pressure", false},
  {OperationKind::calc_temperature, "calc-temperature", false},
  {OperationKind::solve_chemistry, "solve-chemistry-dt=", true},
};

const std::regex flt_regex("^" FLT_PATTERN "$");

/// Convert the command line argument to an OperationSpec
///
/// op_spec_str_view is simply a view of ``arg`` without the leading flag.
OperationSpec op_spec_from_arg_(const char* arg,
                                std::string_view op_spec_str_view) {
  

  for (const OperationProps& op_prop : op_prop_list){
    if (op_prop.requires_dt && starts_with(op_spec_str_view, op_prop.name)){
      std::string dt_str = std::string(
          op_spec_str_view.substr(std::strlen(op_prop.name))
      );
      if (!std::regex_match(dt_str, flt_regex)) {
        std::fprintf(stderr,
            "problem parsing `%s`. Were you trying to specify `-O%s<float>`?\n",
            arg, op_prop.name);
        std::exit(1);
      } 
      return {op_prop.kind, std::stod(dt_str)};
    } else if ((!op_prop.requires_dt) && op_spec_str_view == op_prop.name) {
      return {op_prop.kind, -1.0};
    }
  }

  // if we made it this far, then there is an error
  std::fprintf(
      stderr, "the `%s` argument corresponds to an unknown operations\n",
      arg);
  std::exit(1);
}

} // anonymous namespace


bool try_parse_op_spec(const char* arg, std::optional<OperationSpec>& op_spec)
{
  if (!starts_with(arg, "-O")) return false;
  std::string_view op_str_view(arg+2);
  OperationSpec tmp = op_spec_from_arg_(arg, op_str_view);

  if (op_spec.has_value()) {
    std::fprintf(stderr, "`%s` is the 2nd arg to specify an operation\n",
                 arg);
    std::exit(1);
  }

  std::optional<OperationSpec> rslt{tmp};
  op_spec = rslt;
  return true;
}

std::string string_repr(const OperationSpec& op_spec) {
  for (const OperationProps& op_prop : op_prop_list) {
    if (op_prop.kind == op_spec.kind) {
      if (!op_prop.requires_dt) return op_prop.name;
      return op_prop.name + std::to_string(op_spec.dt);
    }
  }
  GRCLI_ERROR(
      "Something went wrong while converting an operation spec to a string. "
      "Somehow we didn't recognize the OperationKind. Did we add a new one "
      "recently and forget to update op_prop_list?"
  );
}
