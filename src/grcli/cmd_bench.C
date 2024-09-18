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

#ifdef USE_BENCHMARK
#include <benchmark/benchmark.h>
#endif

#include "args/CliParser.h"
#include "args/EscapedArgsConfig.h"
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
  args::BackendConfig backend_arg_spec;

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
    } else if (backend_arg_spec.try_parse(ptr, parser)) {
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
  FieldData fields = scenario::initialize_grid(*pack.get_chemistry_data(),
                                               pack.initial_units(),
                                               grid_spec);

  GrackleDriver driver(pack.get_chemistry_data(), pack.chemistry_storage(),
                       pack.initial_units(), std::move(fields),
                       op_spec.value());

  if (false) {
    int n_iter = 20;
    grackle::BenchState tmp(n_iter);
    driver(tmp);

    double n_seconds = tmp.GetTotalElapsedSeconds();
    std::printf("%d iterations, %g seconds\n", n_iter, n_seconds);
  } else {

    auto my_test = [&driver](benchmark::State& st) { driver(st); };
    benchmark::RegisterBenchmark(string_repr(op_spec.value()), my_test);

    // construct argv and argc
    // ToDo: there is probably a smarter way to do this without copying the
    //       vector (since we aren't copying the pointers this isn't terrible)
    std::vector<char*> argv_vec = backend_arg_spec.args_copy();
    int argc = int(argv_vec.size()) + 1;
    argv_vec.reserve(argc + 2);
    argv_vec.insert(argv_vec.begin(), parser.bin_name());
    argv_vec.push_back(nullptr);

    benchmark::Initialize(&argc, (char**)(argv_vec.data()));
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();

  }
  std::exit(0);
}

