#include <limits>
#include <utility>  // std::move, std::pair

#include <benchmark/benchmark.h>

#include <grackle.h>
#include "grtestutils/GrackleCtxPack.hpp"
#include "grtestutils/field_container.hpp"
#include "grtestutils/fill_field_vals.hpp"
#include "grtestutils/preset.hpp"
#include "grtestutils/status.hpp"

// define the benchmark
static void BM_solve_chemistry(benchmark::State& state) {
  int64_t n_contiguous_cells = state.range(0);
  if (n_contiguous_cells > std::numeric_limits<int>::max()) {
    state.SkipWithError("number of cells along contiguous axis is too big");
    return;  // Early return is allowed when SkipWithError() has been used.
  }

  // Perform setup here
  // setup the context for a grackle solver with metal_cooling
  grtest::ParamConf conf(grtest::ChemPreset::primchem0,
                         grtest::InitialUnitPreset::simple_z0,
                         grtest::make_ParamPair_vec({
                             {"dust_chemistry", 0},
                         }));
  std::pair<grtest::GrackleCtxPack, grtest::Status> ctx_rslt =
      grtest::GrackleCtxPack::create(conf);
  if (ctx_rslt.second.is_err()) {
    state.SkipWithError("Error while creating grtest::GrackleCtxPack: " +
                        ctx_rslt.second.to_string());
    return;  // Early return is allowed when SkipWithError() has been used.
  }
  grtest::GrackleCtxPack ctx_pack = std::move(ctx_rslt.first);

  // we reuse the initial units as the current units
  code_units my_units_obj = ctx_pack.initial_units();

  // make and fill the field container
  int field_grid_dims[3] = {static_cast<int>(n_contiguous_cells), 1, 1};
  std::pair<grtest::GridLayout, grtest::Status> layout_rslt =
      grtest::GridLayout::try_from_dim(1, field_grid_dims);
  if (layout_rslt.second.is_err()) {
    state.SkipWithError("Error while creating grtest::GridLayout: " +
                        ctx_rslt.second.to_string());
    return;  // Early return is allowed when SkipWithError() has been used.
  }
  std::pair<grtest::FieldContainer, grtest::Status> tmp =
      grtest::create_and_fill_FieldContainer(
          grtest::make_simple_tile(ctx_pack, my_units_obj, 1.0), ctx_pack,
          layout_rslt.first);
  if (tmp.second.is_err()) {
    state.SkipWithError(
        "Error while creating/filling grtest::FieldContainer: " +
        ctx_rslt.second.to_string());
    return;  // Early return is allowed when SkipWithError() has been used.
  }
  const grtest::FieldContainer fc_initial = std::move(tmp.first);

  // create the fluid-container with buffers that will be modified during the
  // calculation
  grtest::FieldContainer fc = fc_initial.clone();

  chemistry_data* my_chemistry = ctx_pack.my_chemistry();
  chemistry_data_storage* my_rates = ctx_pack.my_rates();
  code_units* my_units = &my_units_obj;
  grackle_field_data* my_fields = fc.get_ptr();

  const double dt_val = 3.15e7 * 1e6 / my_units->time_units;

  // the following loop is actually benchmarked
  for (auto _ : state) {
    // pause the timer so that we can reset the state of fc
    // -> google-benchmark docs do NOT recomended this (since it incurs a lot
    //    of overhead)
    // -> but it is unavoidable
    state.PauseTiming();
    // reset the values in fc to the desired initial field values
    // (this is an extremely quick operation)
    fc_initial.copy_into(fc, /*bypass_check=*/true);
    state.ResumeTiming();

    // the actual function that is profiled
    local_solve_chemistry(my_chemistry, my_rates, my_units, my_fields, dt_val);
  }
}

// Register the function as a benchmark
BENCHMARK(BM_solve_chemistry)
    ->Name("solve_chemistry")
    ->Arg(1)
    ->Arg(8)
    ->Arg(512);

// Run the benchmark
BENCHMARK_MAIN();
