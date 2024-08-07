/***********************************************************************
/
/ Implement function used to run the bench subcommand
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <string_view>
#include <optional>

#include "CliParser.h"
#include "cmd_bench.h"
#include "executor.h"
#include "grid_problem.h"
#include "handle_gr_parameters.h"
#include "operation.h"
#include "utils.h"

/// execute the bench command
[[noreturn]] void cmd::bench::run(CliParser& parser)
{
  std::optional<CliParamSpec> gr_param_spec;
  scenario::CliGridSpec grid_spec;
  std::optional<OperationSpec> op_spec;

  const char* ptr;
  while ((ptr = parser.next()) != nullptr) {
    std::string_view arg(ptr);
    try_parse_help(ptr, parser.bin_name());
    if (try_parse_cli_paramspec(ptr, parser, gr_param_spec)) {
      continue;
    } else if (try_parse_cli_grid_component(ptr, grid_spec)) {
      continue;
    } else if (try_parse_op_spec(ptr, op_spec)) {
      continue;
    }
    err_unrecognized_arg(ptr);
  }


  // move onto the actual work:

  if (!op_spec.has_value()) {
    std::fprintf(stderr, "The operation wasn't specified");
    std::exit(1);
  }

  // initialize the full grackle solver
  FullGrackleSolverPack pack = create_full_grackle_solver(gr_param_spec);

  // initialize the field data
  if (!grid_spec.any_specified) {
    std::fprintf(stderr, "The scenario wasn't specified");
    std::exit(1);
  }
  FieldData fields = scenario::initialize_grid(*pack.chemistry_data(),
                                               pack.initial_units(),
                                               grid_spec);


  GRCLI_ERROR("NOT FINISHED IMPLEMENTING SUBCOMMAND YET");
  std::exit(1);
}

