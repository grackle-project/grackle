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
#include "grid_problem.h"
#include "handle_gr_parameters.h"
#include "utils.h"

/// execute the bench command
[[noreturn]] void cmd::bench::run(CliParser& parser)
{
  std::optional<CliParamSpec> gr_param_spec;

  const char* ptr;
  while ((ptr = parser.next()) != nullptr) {
    try_parse_help(ptr, parser.bin_name());
    if (try_parse_cli_paramspec(ptr, parser, gr_param_spec)) continue;
    err_unrecognized_arg(ptr);
  }

  FullGrackleSolverPack pack = create_full_grackle_solver(gr_param_spec);

  GRCLI_ERROR("NOT FINISHED IMPLEMENTING SUBCOMMAND YET");
  std::exit(1);
}

