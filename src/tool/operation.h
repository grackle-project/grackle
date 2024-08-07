#ifndef GRCLI_OPERATION_H
#define GRCLI_OPERATION_H

#include <optional>

enum struct OperationKind {
  calc_cooling_time,
  calc_dust_temperature,
  calc_pressure,
  calc_temperature,
  solve_chemistry
};

struct OperationSpec {
  OperationKind kind;
  double dt;
};

/// Checks whether the specified arg correspond to an OperationSpec:
///
/// * If so, this will parse the argument, update the op_spec
///   argument, and this function will return true
///
/// * Otherwise, this returns false
bool try_parse_op_spec(const char* arg, std::optional<OperationSpec>& op_spec);

#endif /* GRCLI_OPERATION_H */
