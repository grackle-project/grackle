#include <cstdio>
#include <cstdlib>
#include <regex>
#include <string_view>
#include <string>

#include "operation.h"
#include "utils.h"

static const std::regex flt_regex("^" FLT_PATTERN "$");

bool try_parse_op_spec(const char* arg, std::optional<OperationSpec>& op_spec)
{
  if (!starts_with(arg, "-O")) return false;

  std::string_view arg_view(arg);

  OperationKind kind;
  double dt = -1.0;
  if (arg_view == "-Ocalc-cooling-time") {
    kind = OperationKind::calc_cooling_time;
  } else if (arg_view == "-Ocalc-dust-temperature") {
    kind = OperationKind::calc_dust_temperature;
  } else if (arg_view == "-Ocalc-pressure") {
    kind = OperationKind::calc_pressure;
  } else if (arg_view == "-Ocalc-temperature") {
    kind = OperationKind::calc_temperature;
  } else if (starts_with(arg_view, "-Osolve-chemistry-dt=")) {
    kind = OperationKind::solve_chemistry;
    std::string dt_str = std::string(arg_view.substr(21));
    if (!std::regex_match(dt_str, flt_regex)) {
      std::fprintf(stderr,
          "problem parsing `%s`. Were you trying to specify "
          "`-Osolve-chemistry-dt=<float>`?\n", arg);
      std::exit(1);
    }
    dt = std::stod(dt_str);
  } else {
    std::fprintf(
        stderr, "the `%s` argument corresponds to an unknown operations\n",
        arg);
    std::exit(1);
  }

  if (op_spec.has_value()) {
    std::fprintf(stderr, "`%s` is the 2nd arg to specify an operation\n",
                 arg);
    std::exit(1);
  }

  std::optional<OperationSpec> rslt{OperationSpec{kind,dt}};
  op_spec = rslt;
  return true;
}
