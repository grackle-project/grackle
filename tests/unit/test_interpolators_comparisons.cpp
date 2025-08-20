#include <random>   // std::minstd_rand
#include <cstdint>  // std::int32_t
#include <cstdio>   // std::printf

#include <math.h>

#include "utest_helpers.hpp"
#include "interpolate.hpp"
#include "fortran_func_wrappers.hpp"

// this file introduces tests that directly compare the newer versions of
// interpolation routines against the original Fortran versions
// - this was originally written without using gtest. The tests have all been
//   retrofitted so that they work with gtest
// - consequentially, they may not provide the best example for what a unit test
//   should look like.

// NOTE: this function has been backported to help with unit-testing.
// Eventually, it will be implemented to help define the C version of
// cool1d_cloudy_g
//
// helper function that retrieves index for redshift dimension (of cloudy
// tables) via bisection
// - the index is one-indexed
// - the names of variables have not been changed for backwards compatibility
//   (it may seem counter-intuitive that clGridDim[1] gives the length of
//    clPar2, but that's because in Fortran you would access clGridDim(2) )
// - NOTE: since we define this function in a header, we must declare it as
//   static inline (in C++ we could just declare it as inline)
static inline long long find_zindex(double zr, long long clGridRank,
                                    const long long* clGridDim,
                                    const double* clPar2) {
  if (clGridRank > 2) {
    long long zindex;
    if (zr <= clPar2[0]) {
      zindex = 1;
    } else if (zr >= clPar2[clGridDim[1] - 2]) {
      zindex = clGridDim[1];
    } else if (zr >= clPar2[clGridDim[1] - 3]) {
      zindex = clGridDim[1] - 2;
    } else {
      zindex = 1;
      long long zhighpt = clGridDim[1] - 2;
      while ((zhighpt - zindex) > 1) {
        long long zmidpt = (long long)((zhighpt + zindex) / 2);
        if (zr >= clPar2[zmidpt - 1]) {
          zindex = zmidpt;
        } else {
          zhighpt = zmidpt;
        }
      }
    }
    return zindex;
  } else {
    return 1;
  }
}

/// The idea here is to encapsulate an DataTable that can be used to execute
/// one of the interpolation functions
class InterpTable {
public:
  // delete default constructor
  InterpTable() = delete;

  InterpTable(std::vector<std::vector<double>> param_vals,
              std::vector<double> dataField) {
    if ((param_vals.size() == 0) || (param_vals.size() > 5)) {
      error("InterpTable::InterpTable()",
            "param_vals must have 1 to 5 entries");
    }

    std::vector<gr_i64> gridDim;
    for (std::size_t i = 0; i < param_vals.size(); i++) {
      std::size_t cur_dim_size = param_vals[i].size();
      if (cur_dim_size < 2) {
        error("InterpTable::InterpTable()",
              "each parameter must have 2+ entries");
      }
      gridDim.push_back(cur_dim_size);
    }

    this->gridDim_ = gridDim;
    this->param_vals_ = param_vals;
    this->dataField_ = dataField;

    // sanity check:
    if (this->dataSize() != (gr_i64)dataField_.size()) {
      error("InterpTable::InterpTable()",
            "dataField has the wrong number of entries.");
    }
  }

  const gr_i64* gridDim() const { return gridDim_.data(); }

  const double* gridField() const { return dataField_.data(); }

  gr_i64 dataSize() const {
    gr_i64 out = 1;
    for (auto elem : this->gridDim_) {
      out *= elem;
    }
    return out;
  }

  int ndim() const { return static_cast<int>(gridDim_.size()); }

  const double* gridPar_0Indexed(int dim) const {
    if (dim < 0) {
      error("InterpTable::gridPar_0Indexed", "can't have negative dimension");
    } else if (dim >= ndim()) {
      error("InterpTable::gridPar_0Indexed", "dim is too large");
    }
    return param_vals_[dim].data();
  }

  double dgridPar_0Indexed(int dim) const {
    const double* gridPar = gridPar_0Indexed(dim);
    gr_i64 gridPar_len = this->gridDim_[dim];

    return ((gridPar[gridPar_len - 1] - gridPar[0]) /
            static_cast<double>(gridPar_len - 1));
  }

private:
  // these shouldn't be mutated after construction
  std::vector<gr_i64> gridDim_;
  // in principle param_vals_ should monotonically increase (and in almost all
  // cases the spacing between values should be constant)
  std::vector<std::vector<double>> param_vals_;
  std::vector<double> dataField_;
};

Timer run_interp_helper(const double* val_vec, const InterpTable& table,
                        std::size_t length, bool use_fortran, double* out) {
  const int rank = table.ndim();

  const double* gridPar1 = table.gridPar_0Indexed(0);
  const double* gridPar2 = (rank > 1) ? table.gridPar_0Indexed(1) : nullptr;
  const double* gridPar3 = (rank > 2) ? table.gridPar_0Indexed(2) : nullptr;
  const double* gridPar4 = (rank > 3) ? table.gridPar_0Indexed(3) : nullptr;
  const double* gridPar5 = (rank > 4) ? table.gridPar_0Indexed(4) : nullptr;

  const double dgridPar1 = table.dgridPar_0Indexed(0);
  const double dgridPar2 = (rank > 1) ? table.dgridPar_0Indexed(1) : 0.0;
  const double dgridPar3 = (rank > 2) ? table.dgridPar_0Indexed(2) : 0.0;
  const double dgridPar4 = (rank > 3) ? table.dgridPar_0Indexed(3) : 0.0;
  const double dgridPar5 = (rank > 4) ? table.dgridPar_0Indexed(4) : 0.0;

  const gr_i64* gridDim = table.gridDim();
  const gr_i64 dataSize = table.dataSize();
  const double* dataField = table.gridField();

  Timer t;

  if ((rank == 1) && (use_fortran)) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
          val_vec[i], gridDim, gridPar1, dgridPar1, dataSize, dataField);
    }
    t.stop();

  } else if (rank == 1) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_1d_g(val_vec[i], gridDim, gridPar1,
                                               dgridPar1, dataSize, dataField);
    }
    t.stop();

  } else if ((rank == 2) && use_fortran) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
          val_vec[i * 2], val_vec[i * 2 + 1], gridDim, gridPar1, dgridPar1,
          gridPar2, dgridPar2, dataSize, dataField);
    }
    t.stop();

  } else if ((rank == 2)) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_2d_g(
          val_vec[i * 2], val_vec[i * 2 + 1], gridDim, gridPar1, dgridPar1,
          gridPar2, dgridPar2, dataSize, dataField);
    }
    t.stop();

  } else if ((rank == 3) && use_fortran) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_3d_g(
          val_vec[i * 3], val_vec[i * 3 + 1], val_vec[i * 3 + 2], gridDim,
          gridPar1, dgridPar1, gridPar2, dgridPar2, gridPar3, dgridPar3,
          dataSize, dataField);
    }
    t.stop();

  } else if (rank == 3) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_3d_g(
          val_vec[i * 3], val_vec[i * 3 + 1], val_vec[i * 3 + 2], gridDim,
          gridPar1, dgridPar1, gridPar2, dgridPar2, gridPar3, dgridPar3,
          dataSize, dataField);
    }
    t.stop();

  } else if ((rank == 4) && use_fortran) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_4d_g(
          val_vec[i * 4], val_vec[i * 4 + 1], val_vec[i * 4 + 2],
          val_vec[i * 4 + 3], gridDim, gridPar1, dgridPar1, gridPar2, dgridPar2,
          gridPar3, dgridPar3, gridPar4, dgridPar4, dataSize, dataField);
    }
    t.stop();

  } else if (rank == 4) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_4d_g(
          val_vec[i * 4], val_vec[i * 4 + 1], val_vec[i * 4 + 2],
          val_vec[i * 4 + 3], gridDim, gridPar1, dgridPar1, gridPar2, dgridPar2,
          gridPar3, dgridPar3, gridPar4, dgridPar4, dataSize, dataField);
    }
    t.stop();

  } else if (use_fortran) {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_5d_g(
          val_vec[i * 5], val_vec[i * 5 + 1], val_vec[i * 5 + 2],
          val_vec[i * 5 + 3], val_vec[i * 5 + 4], gridDim, gridPar1, dgridPar1,
          gridPar2, dgridPar2, gridPar3, dgridPar3, gridPar4, dgridPar4,
          gridPar5, dgridPar5, dataSize, dataField);
    }
    t.stop();

  } else {
    t.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_5d_g(
          val_vec[i * 5], val_vec[i * 5 + 1], val_vec[i * 5 + 2],
          val_vec[i * 5 + 3], val_vec[i * 5 + 4], gridDim, gridPar1, dgridPar1,
          gridPar2, dgridPar2, gridPar3, dgridPar3, gridPar4, dgridPar4,
          gridPar5, dgridPar5, dataSize, dataField);
    }
    t.stop();
  }
  return t;
}

struct interp_3dz_arg_ {
  double par1;
  double zr;
  double par3;
  gr_i64 zindex;
  gr_i64 end_int;
};

Timer run_interp_3dz_(const double* vals, const InterpTable& table,
                      std::size_t length, bool use_fortran, double* out) {
  const int rank = table.ndim();
  if (rank != 3) error("run_interp_3dz", "requires a rank 3 table");

  const double* gridPar1 = table.gridPar_0Indexed(0);
  const double* gridPar2 = table.gridPar_0Indexed(1);
  const double* gridPar3 = table.gridPar_0Indexed(2);

  const double dgridPar1 = table.dgridPar_0Indexed(0);
  const double dgridPar3 = table.dgridPar_0Indexed(2);

  const gr_i64* gridDim = table.gridDim();
  const gr_i64 dataSize = table.dataSize();
  const double* dataField = table.gridField();

  // need to do some conversions on the arguments (we try to do this ahead of
  // time to minimize impact on the timings). Ideally, we'd do this before
  // calling this function so that we could reuse the same array for both
  // versions of the function
  std::vector<interp_3dz_arg_> inputs_vec(length);
  for (std::size_t i = 0; i < length; i++) {
    const double zr = vals[i * 3 + 1];  // z redshift

    // Calculate index for redshift dimension - intentionally kept 1-indexed
    const long long zindex = find_zindex(zr, rank, gridDim, gridPar2);
    const gr_i64 end_int = ((rank > 2) && (zindex == gridDim[1]));

    inputs_vec[i] = {vals[i * 3], zr, vals[i * 3 + 2], zindex, end_int};
  }

  const interp_3dz_arg_* inputs = inputs_vec.data();

  Timer timer;

  if (use_fortran) {
    timer.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::fortran_wrapper::interpolate_3dz_g(
          inputs[i].par1, inputs[i].zr, inputs[i].par3, gridDim, gridPar1,
          dgridPar1, gridPar2, (inputs[i].zindex), gridPar3, dgridPar3,
          dataSize, dataField, inputs[i].end_int);
    }
    timer.stop();

  } else {
    timer.start();
    for (std::size_t i = 0; i < length; i++) {
      out[i] = grackle::impl::interpolate_3dz_g(
          inputs[i].par1, inputs[i].zr, inputs[i].par3, gridDim, gridPar1,
          dgridPar1, gridPar2, inputs[i].zindex, gridPar3, dgridPar3, dataSize,
          dataField, inputs[i].end_int);
    }
    timer.stop();
  }

  return timer;
}

std::vector<double> run_interp(const std::vector<std::vector<double>>& vals,
                               const InterpTable& table, bool use_fortran,
                               bool use_3dz, double* elapsed = nullptr) {
  const std::size_t ndim = table.ndim();
  std::vector<double> flat_vec(vals.size() * ndim);
  for (std::size_t i = 0; i < vals.size(); i++) {
    for (std::size_t j = 0; j < ndim; j++) {
      flat_vec[i * ndim + j] = vals[i][j];
    }
  }

  std::vector<double> out(vals.size());
  Timer t;
  if (use_3dz) {
    t = run_interp_3dz_(flat_vec.data(), table, vals.size(), use_fortran,
                        out.data());
  } else {
    t = run_interp_helper(flat_vec.data(), table, vals.size(), use_fortran,
                          out.data());
  }

  if (elapsed != nullptr) *elapsed = t.get();
  return out;
}

void compare_interp_impls_(const std::vector<std::vector<double>>& inputs,
                           const InterpTable& table, bool use_3dz,
                           bool show_timing = false) {
  std::vector<double> ref, actual;
  double fortran_duration_sec, c_duration_sec;

  ref = run_interp(inputs, table, true, use_3dz);  // run an extra-time for
                                                   // more-fair timing
  ref = run_interp(inputs, table, true, use_3dz, &fortran_duration_sec);
  actual = run_interp(inputs, table, false, use_3dz, &c_duration_sec);
  if (show_timing) {
    printf("      ...Timing Compare: Fortran = %g sec C = %g sec\n",
           fortran_duration_sec, c_duration_sec);
  }

  if (false) {  // debugging statement:
    std::string ref_s = vec_to_string(ref);
    std::string actual_s = vec_to_string(actual);
    std::printf("   reference: %s\n", ref_s.c_str());
    std::printf("   actual: %s\n", actual_s.c_str());
  }
  compare_values(actual, ref, 0.0, 0.0, "");
};

// Miscelaneous Functions
// ======================

// uniform distribution on the interval (0, 1]
double uniform_dist_transform_(std::minstd_rand& prng) {
  // sanity check to confirm that the largest value returned by generator is
  // representable (without any loss of precision) by a double (its <= 2^53)
  static_assert((std::minstd_rand::max() <= 9007199254740992) &
                    (std::minstd_rand::min() == 1),
                "sanity check failed");

  return (static_cast<double>(prng()) / static_cast<double>(prng.max()));
}

// Functions for constructing InterpTable
// ======================================

struct ParamProps {
  double min_val;
  double max_val;
  int num_vals;
};

InterpTable build_table(std::uint32_t seed,
                        const std::vector<ParamProps>& paramprop_vec) {
  std::minstd_rand prng = std::minstd_rand(seed);
  std::vector<std::vector<double>> param_vals;
  std::size_t field_size = 1;

  for (const ParamProps& cur_param_prop : paramprop_vec) {
    const int dim_size = cur_param_prop.num_vals;

    // initialize cur_param_vals
    std::vector<double> cur_param_vals(dim_size);
    double dt = (cur_param_prop.max_val - cur_param_prop.min_val) / dim_size;
    for (int i = 0; i < dim_size; i++) {
      cur_param_vals[i] = cur_param_prop.min_val + i * dt;
    }
    cur_param_vals[dim_size - 1] = cur_param_prop.max_val;
    param_vals.push_back(cur_param_vals);

    // update field_size
    field_size *= dim_size;
  }

  // now initialize field_data
  std::vector<double> field_data(field_size);
  for (std::size_t i = 0; i < field_size; i++) {
    field_data[i] = 100.0 * uniform_dist_transform_(prng) - 50.0;
  }

  return InterpTable(param_vals, field_data);
}

InterpTable get_3dz_table() {
  DummyGrackleConfig config(3, 0.0, true, true);

  cloudy_data& cloud_data_obj = config.chem_rates.cloudy_primordial;

  const int rank = (int)cloud_data_obj.grid_rank;

  // setup param_vals
  std::vector<std::vector<double>> param_vals;
  for (int i = 0; i < rank; i++) {
    const int len = (int)(cloud_data_obj.grid_dimension[i]);
    std::vector<double> cur(len);
    for (int j = 0; j < len; j++) {
      cur[j] = cloud_data_obj.grid_parameters[i][j];
    }
    param_vals.push_back(cur);
  }

  // setup dataField
  const long long data_size = (long long)cloud_data_obj.data_size;
  std::vector<double> dataField(data_size);
  for (long long i = 0; i < data_size; i++) {
    dataField[i] = cloud_data_obj.heating_data[i];
  }

  return {param_vals, dataField};
}

// Functions to help generate locations to perform interpolations at
// =================================================================

// this returns a vector of vectors appropriate interpolate function for cases
// where 1 or more inputs come from special_vals (the other inputs are taken
// from ordinary_vals)
//
// EXAMPLES
// --------
// if special_vals = {1.0,2.0} & ordinary_vals = {-1.0, -2.0}, this passes the
// following to run_interp in separate calls:
//     { 1.0, 2.0}, {-1.0, 2.0}, { 1.0,-2.0}
// if special_vals = {1.0,2.0,3.0} & ordinary_vals = {-1.0, -2.0, -3.0}, this
// passes the following to run_interp in separate calls:
//     { 1.0, 2.0, 3.0}, {-1.0, 2.0, 3.0}, { 1.0,-2.0, 3.0}, {-1.0,-2.0, 3.0}
//     { 1.0, 2.0,-3.0}, {-1.0, 2.0,-3.0}, { 1.0,-2.0,-3.0}
//
// The equivalent in python would be the result of the following command:
// >>> list(itertools.product(*zip(special_vals, ordinary_vals)))[:-1]
std::vector<std::vector<double>> get_combinations_(
    const std::vector<double>& special_vals,
    const std::vector<double>& ordinary_vals) {
  std::size_t rank = special_vals.size();

  int choice[5] = {0, 0, 0, 0, 0};

  std::size_t out_count = std::size_t(pow(2, rank) - 1);
  std::vector<std::vector<double>> out;
  out.reserve(out_count - 1);

  while (out.size() < out_count) {
    std::vector<double> cur(rank);
    for (std::size_t i = 0; i < rank; i++) {
      cur[i] = (choice[i] == 0) ? special_vals[i] : ordinary_vals[i];
    }
    out.push_back(cur);

    // now increment choice
    for (std::size_t i = 0; i < rank; i++) {
      if (choice[i] == 0) {
        choice[i] = 1;
        break;
      }
      choice[i] = 0;
    }
  }

  return out;
}

struct TableSummaryProps {
  std::vector<double> min_param_val;
  std::vector<double> max_param_val;
  std::vector<double> unaligned_inrange;
  std::vector<double> unaligned_inrange_alt;
  std::vector<double> half_min_param_val;
  std::vector<double> double_max_param_val;
};

TableSummaryProps get_interp_table_summary(const InterpTable& table) {
  const int rank = table.ndim();

  // get useful values for test:
  std::vector<double> min_param_val(rank);
  std::vector<double> max_param_val(rank);
  std::vector<double> unaligned_inrange(rank);
  std::vector<double> unaligned_inrange_alt(rank);
  std::vector<double> half_min_param_val(rank);
  std::vector<double> double_max_param_val(rank);
  for (int i = 0; i < rank; i++) {
    const double* cur_param_vals = table.gridPar_0Indexed(i);
    const int cur_param_len = table.gridDim()[i];

    min_param_val[i] = cur_param_vals[0];
    max_param_val[i] = cur_param_vals[cur_param_len - 1];

    double delta = cur_param_vals[1] - cur_param_vals[0];

    unaligned_inrange[i] = cur_param_vals[0] + delta * 0.3;
    unaligned_inrange_alt[i] = cur_param_vals[0] + delta * 0.5;

    half_min_param_val[i] = min_param_val[i] / 2;
    double_max_param_val[i] = max_param_val[i] * 2;
  }

  TableSummaryProps out;
  out.min_param_val = min_param_val;
  out.max_param_val = max_param_val;
  out.unaligned_inrange = unaligned_inrange;
  out.unaligned_inrange_alt = unaligned_inrange_alt;
  out.half_min_param_val = half_min_param_val;
  out.double_max_param_val = double_max_param_val;
  return out;
}

// Define the actual tests
// =======================

void run_test(const InterpTable& table, const bool use_3dz = false) {
  TableSummaryProps props = get_interp_table_summary(table);

  // define the actual tests in the test suite

  std::printf(" -> comparing some in-range values\n");
  {
    std::vector<std::vector<double>> inputs =
        get_combinations_(props.unaligned_inrange, props.unaligned_inrange_alt);
    compare_interp_impls_(inputs, table, use_3dz);
  }

  std::printf(" -> comparing cases with 1+ inputs on left grid edge\n");
  {
    std::vector<std::vector<double>> inputs =
        get_combinations_(props.min_param_val, props.unaligned_inrange);
    compare_interp_impls_(inputs, table, use_3dz);
  }

  std::printf(" -> comparing cases with 1+ inputs on right grid edge\n");
  {
    std::vector<std::vector<double>> inputs =
        get_combinations_(props.max_param_val, props.unaligned_inrange);
    compare_interp_impls_(inputs, table, use_3dz);
  }

  if (!use_3dz) {
    std::printf(" -> comparing cases with 1+ inputs less than min grid val\n");
    std::vector<std::vector<double>> inputs =
        get_combinations_(props.half_min_param_val, props.unaligned_inrange);
    compare_interp_impls_(inputs, table, use_3dz);
  }

  std::printf(" -> comparing cases with 1+ inputs greater than max grid val\n");
  {
    std::vector<std::vector<double>> inputs =
        get_combinations_(props.double_max_param_val, props.unaligned_inrange);
    compare_interp_impls_(inputs, table, use_3dz);
  }

  if (use_3dz) {
    std::printf(" -> comparing cases with some extra redshifts\n");

    gr_i64 n_gridded_z_vals = table.gridDim()[1];

    double temp[5] = {table.gridPar_0Indexed(1)[n_gridded_z_vals - 4],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 3],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 2],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 1],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 1] * 1.01};

    std::vector<std::vector<double>> inputs;
    for (int i = 0; i < 5; i++) {
      inputs.push_back(
          {props.unaligned_inrange[0], temp[i], props.unaligned_inrange[2]});
    }
    compare_interp_impls_(inputs, table, use_3dz);
  }

  std::size_t n_random = 1000000;
  std::printf(" -> comparing %d random points within grid (with timing)\n",
              int(n_random));
  {
    const int rank = table.ndim();
    // the fact that we use vectors of vectors as inputs makes this inefficient
    std::minstd_rand prng = std::minstd_rand(353545);
    std::vector<std::vector<double>> inputs(n_random);
    for (std::size_t i = 0; i < n_random; i++) {
      std::vector<double> tmp(rank);
      for (int j = 0; j < rank; j++) {
        double max_val = props.max_param_val[j];
        double min_val = props.min_param_val[j];
        tmp[j] = uniform_dist_transform_(prng) * (max_val - min_val) + min_val;
      }
      inputs[i] = tmp;
    }

    compare_interp_impls_(inputs, table, use_3dz, true);
  }
}

const std::uint32_t seed = 342;

TEST(InterpolationTest, CompareInterpolate1D) {
  std::printf("comparing interpolate_1d_g:\n");
  InterpTable table = build_table(seed, {{-6.0, 6.0, 25}});
  run_test(table);
}

TEST(InterpolationTest, CompareInterpolate2D) {
  std::printf("\ncomparing interpolate_2d_g:\n");
  InterpTable table = build_table(seed, {{-6.0, 6.0, 25}, {0.0, 10.0, 11}});
  run_test(table);
}

TEST(InterpolationTest, CompareInterpolate3D) {
  std::printf("\ncomparing interpolate_3d_g:\n");
  InterpTable table =
      build_table(seed, {{-6.0, 6.0, 25}, {0.0, 10.0, 11}, {-1.0, 0.0, 5}});
  run_test(table);
}

TEST(InterpolationTest, CompareInterpolate4D) {
  std::printf("\ncomparing interpolate_4d_g:\n");
  InterpTable table = build_table(
      seed, {{-6.0, 6.0, 25}, {0.0, 10.0, 11}, {-1.0, 0.0, 5}, {0.5, 5.5, 6}});
  run_test(table);
}

TEST(InterpolationTest, CompareInterpolate5D) {
  std::printf("\ncomparing interpolate_5d_g:\n");
  InterpTable table = build_table(seed, {{-6.0, 6.0, 25},
                                         {0.0, 10.0, 11},
                                         {-1.0, 0.0, 5},
                                         {0.5, 5.5, 6},
                                         {-10.0, 0.0, 11}});
  run_test(table);
}

TEST(InterpolationTest, CompareInterpolate3Dz) {
  std::printf("\ncomparing interpolate_3dz_g:\n");
  InterpTable table = get_3dz_table();
  run_test(table, true);
}
