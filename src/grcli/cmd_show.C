/***********************************************************************
/
/ Declare functions used to run the show-* subcommands.
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <cstdlib>
#include <cstdio>
#include <optional>

#include "../clib/show.h"

#include "cmd_show.h"
#include "handle_gr_parameters.h"
#include "utils.h"



[[noreturn]] void cmd::show::run(CliParser& parser, cmd::show::Kind kind) {

  if (kind == cmd::show::Kind::initial_units) {

    if (parser.has_next()) err_unrecognized_arg(parser.next());
    code_units initial_units = get_default_units();
    show_units(stdout, &initial_units);

  } else if (kind == cmd::show::Kind::parameters) {

    std::optional<CliParamSpec> gr_param_spec;
    const char* ptr;
    while ((ptr = parser.next()) != nullptr) {
      try_parse_help(ptr, parser.bin_name());
      if (try_parse_cli_paramspec(ptr, parser, gr_param_spec)) continue;
      err_unrecognized_arg(ptr);
    }
    FullGrackleSolverPack tmp = create_full_grackle_solver(gr_param_spec);
    show_parameters(stdout,  tmp.chemistry_data());

  } else if (kind == cmd::show::Kind::precision) {

    if (parser.has_next()) err_unrecognized_arg(parser.next());
    std::printf("float-point-size = %d bytes\n", int(sizeof(gr_float)));

  } else if (kind == cmd::show::Kind::version) {

    if (parser.has_next()) err_unrecognized_arg(parser.next());
    show_version(stdout);

  } else {
    GRCLI_ERROR("We forgot to update this function when we added a new enum");
  }

  std::printf("Warning, the format/contents of this command may change\n");
  std::exit(0);
}
