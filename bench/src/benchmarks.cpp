#include <cstdint>
#include <limits>
#include <optional>
#include <utility>  // std::move, std::pair

#include <benchmark/benchmark.h>

#include <grackle.h>
#include "grtestutils/GrackleCtxPack.hpp"
#include "grtestutils/field_container.hpp"
#include "grtestutils/fill_field_vals.hpp"
#include "grtestutils/preset.hpp"
#include "grtestutils/status.hpp"

enum struct FieldIC { UNIFORM };

struct ArgPack {
  FieldIC field_ic;  // <- this is mostly just a placeholder
  grtest::ParamConf param_conf;
  std::optional<int> total_cell_count;
};

struct SetupRslt {
  grtest::GrackleCtxPack ctx_pack;
  code_units unit_obj;
  grtest::FieldContainer fc;
};

static std::pair<std::optional<SetupRslt>, grtest::Status> perform_setup(
    const ArgPack& arg_pack, int64_t n_contiguous_cells) {
  // create a useful shorthand
  using OutPair = std::pair<std::optional<SetupRslt>, grtest::Status>;
  auto fail_obj = [](const std::string& message) -> OutPair {
    return OutPair{std::nullopt, grtest::error::Adhoc(message)};
  };

  // setup the context for a grackle solver
  std::pair<grtest::GrackleCtxPack, grtest::Status> ctx_rslt =
      grtest::GrackleCtxPack::create(arg_pack.param_conf);
  if (ctx_rslt.second.is_err()) {
    return fail_obj("Error while creating grtest::GrackleCtxPack: " +
                    ctx_rslt.second.to_string());
  }
  grtest::GrackleCtxPack ctx_pack = std::move(ctx_rslt.first);

  // we reuse the initial units as the current units
  code_units my_units_obj = ctx_pack.initial_units();

  // make and fill the field container
  if (n_contiguous_cells > std::numeric_limits<int>::max()) {
    return fail_obj("n_contiguous_cells is too big");
  } else if (n_contiguous_cells <= 0) {
    return fail_obj("n_contiguous_cells must exceed 0");
  }
  int n_contig_cells = static_cast<int>(n_contiguous_cells);
  int total_cell_count = arg_pack.total_cell_count.value_or(n_contig_cells);
  int ratio = total_cell_count / n_contig_cells;
  int remainder = total_cell_count % n_contig_cells;
  if (ratio <= 0 || remainder != 0) {
    return fail_obj("n_contiguous_cells must not exceed total_cell_count");
  }

  int field_grid_dims[3] = {n_contig_cells, ratio, 1};
  std::pair<grtest::GridLayout, grtest::Status> layout_rslt =
      grtest::GridLayout::try_from_dim(3, field_grid_dims);
  if (layout_rslt.second.is_err()) {
    return fail_obj("Error while creating grtest::GridLayout: " +
                    ctx_rslt.second.to_string());
  }
  if (arg_pack.field_ic != FieldIC::UNIFORM) {
    return fail_obj("can't handle field_ic argument (yet!)");
  }
  std::pair<grtest::FieldContainer, grtest::Status> tmp =
      grtest::create_and_fill_FieldContainer(
          grtest::make_simple_tile(ctx_pack, my_units_obj, 1.0), ctx_pack,
          layout_rslt.first);
  if (tmp.second.is_err()) {
    return fail_obj("Error while creating/filling grtest::FieldContainer: " +
                    ctx_rslt.second.to_string());
  }

  return OutPair{std::optional<SetupRslt>{
                     {std::move(ctx_pack), my_units_obj, std::move(tmp.first)}},
                 grtest::OkStatus()};
}

// define the benchmark
static void BM_solve_chemistry(benchmark::State& state, const ArgPack pack) {
  int64_t n_contiguous_cells = state.range(0);

  std::pair<std::optional<SetupRslt>, grtest::Status> setup_pair =
      perform_setup(pack, n_contiguous_cells);
  if (setup_pair.second.is_err()) {
    state.SkipWithError("number of cells along contiguous axis is too big");
    return;  // Early return is allowed when SkipWithError() has been used.
  } else if (!setup_pair.first.has_value()) {
    state.SkipWithError("unexpected error!");
    return;  // Early return is allowed when SkipWithError() has been used.
  }

  grtest::GrackleCtxPack ctx_pack = std::move(setup_pair.first->ctx_pack);
  code_units my_units_obj = setup_pair.first->unit_obj;
  const grtest::FieldContainer fc_initial = std::move(setup_pair.first->fc);

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

static ArgPack make_uniform_arg_pack(
    grtest::ChemPreset chem_preset,
    std::optional<int> total_cell_count = std::nullopt) {
  return {FieldIC::UNIFORM,
          grtest::ParamConf(chem_preset, grtest::InitialUnitPreset::simple_z0,
                            grtest::make_ParamPair_vec({
                                {"dust_chemistry", 0},
                            })),
          total_cell_count};
}

// Register the function as a benchmark
static const ArgPack uniform_pc0_pack =
    make_uniform_arg_pack(grtest::ChemPreset::primchem0);
BENCHMARK_CAPTURE(BM_solve_chemistry, UniformPC0, uniform_pc0_pack)
    ->Arg(8)
    ->Arg(64)
    ->Arg(512);

static const ArgPack uniform_pc0_pack_512cells = make_uniform_arg_pack(
    grtest::ChemPreset::primchem0, std::optional<int64_t>{512});
BENCHMARK_CAPTURE(BM_solve_chemistry, 512CellUniformPC0,
                  uniform_pc0_pack_512cells)
    ->Arg(8)
    ->Arg(64);

static const ArgPack uniform_pc1_pack =
    make_uniform_arg_pack(grtest::ChemPreset::primchem1);
BENCHMARK_CAPTURE(BM_solve_chemistry, UniformPC1, uniform_pc1_pack)
    ->Arg(8)
    ->Arg(64)
    ->Arg(512);

static const ArgPack uniform_pc1_pack_512cells = make_uniform_arg_pack(
    grtest::ChemPreset::primchem1, std::optional<int64_t>{512});
BENCHMARK_CAPTURE(BM_solve_chemistry, 512CellUniformPC1,
                  uniform_pc1_pack_512cells)
    ->Arg(8)
    ->Arg(64);

static const ArgPack uniform_pc2_pack =
    make_uniform_arg_pack(grtest::ChemPreset::primchem2);
BENCHMARK_CAPTURE(BM_solve_chemistry, UniformPC2, uniform_pc2_pack)->Arg(64);

static const ArgPack uniform_pc3_pack =
    make_uniform_arg_pack(grtest::ChemPreset::primchem3);
BENCHMARK_CAPTURE(BM_solve_chemistry, UniformPC3, uniform_pc3_pack)->Arg(64);

// Run the benchmark
BENCHMARK_MAIN();
