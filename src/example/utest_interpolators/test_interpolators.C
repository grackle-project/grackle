#include <random> // std::minstd_rand
#include <cstdint> // std::int32_t
#include <cstdio> // std::printf

#include <math.h>

#include "../utest_helpers.hpp"


typedef long long gr_int64;

extern "C" {

#include "../../clib/interop/interop_funcs.h"

void FORTRAN_NAME(interpolate_1d_g)(
        const double* input1, const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_2d_g)(
        const double* input1, const double* input2,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const double* dgridPar2,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_3dz_g)(
        const double* input1, const double* input2, const double* input3,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const gr_int64* index2,
        const double* gridPar3, const double* dgridPar3,
        const gr_int64* dataSize, const double* dataField,
        const gr_int64* end_int,
        double* value);

void FORTRAN_NAME(interpolate_3d_g)(
        const double* input1, const double* input2, const double* input3,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const double* dgridPar2,
        const double* gridPar3, const double* dgridPar3,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_4d_g)(
        const double* input1, const double* input2, const double* input3,
        const double* input4,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const double* dgridPar2,
        const double* gridPar3, const double* dgridPar3,
        const double* gridPar4, const double* dgridPar4,
        const gr_int64* dataSize, const double* dataField,
        double* value);

void FORTRAN_NAME(interpolate_5d_g)(
        const double* input1, const double* input2, const double* input3,
        const double* input4, const double* input5,
        const gr_int64* gridDim,
        const double* gridPar1, const double* dgridPar1,
        const double* gridPar2, const double* dgridPar2,
        const double* gridPar3, const double* dgridPar3,
        const double* gridPar4, const double* dgridPar4,
        const double* gridPar5, const double* dgridPar5,
        const gr_int64* dataSize, const double* dataField,
        double* value);

}



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
                                    const double* clPar2){
  if (clGridRank > 2){
    long long zindex;
    if (zr <= clPar2[0]) {
      zindex = 1;
    } else if (zr >= clPar2[clGridDim[1]-2]) {
      zindex = clGridDim[1];
    } else if (zr >= clPar2[clGridDim[1]-3]) {
      zindex = clGridDim[1] - 2;
    } else {
      zindex = 1;
      long long zhighpt = clGridDim[1] - 2;
      while ((zhighpt - zindex) > 1) {
        long long zmidpt = (long long)((zhighpt + zindex) / 2);
        if (zr >= clPar2[zmidpt-1]){
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



class InterpTable {

public:
  // delete default constructor
  InterpTable() = delete;

  InterpTable(std::vector<std::vector<double>> param_vals,
              std::vector<double> dataField)
  {
    if ((param_vals.size() == 0) || (param_vals.size() > 5)){
      error("InterpTable::InterpTable()",
            "param_vals must have 1 to 5 entries");
    }

    std::vector<gr_int64> gridDim;
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
    if (this->dataSize() != dataField_.size()) {
      error("InterpTable::InterpTable()",
            "dataField has the wrong number of entries.");
    }
  }

  const gr_int64* gridDim() const
  { return gridDim_.data(); }

  const double* gridField() const
  { return dataField_.data(); }

  gr_int64 dataSize() const {
    gr_int64 out = 1;
    for (auto elem: this->gridDim_) { out *= elem; }
    return out;
  }

  int ndim() const
  { return static_cast<int>(gridDim_.size()); }

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
    gr_int64 gridPar_len = this->gridDim_[dim];

    return ((gridPar[gridPar_len - 1] - gridPar[0]) /
            static_cast<double>(gridPar_len - 1));
  }

private:
  // these shouldn't be mutated after construction
  std::vector<gr_int64> gridDim_;
  // in principle param_vals_ should monotonically increase (and in almost all
  // cases the spacing between values should be constant)
  std::vector<std::vector<double>> param_vals_;
  std::vector<double> dataField_;
};

double run_interp(const std::vector<double>& val_vec,
                  const InterpTable& table,
                  bool use_fortran){

  const int rank = table.ndim();
  if (val_vec.size() != rank){
    error("run_interp", "val_vec has wrong number of dims.");
  } else if ((rank <= 0) || (rank > 5)) {
    error("run_interp", "rank has an invalid value.");
  }

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

  const gr_int64* gridDim = table.gridDim();
  const gr_int64 dataSize = table.dataSize();
  const double* dataField = table.gridField();
  
  double out;
  if (rank == 1) {
    if (use_fortran){
      FORTRAN_NAME(interpolate_1d_g)(&val_vec[0], gridDim,
                                     gridPar1, &dgridPar1,
                                     &dataSize, dataField,
                                     &out);
    } else {
      interpolate_1d_g(val_vec[0], gridDim,
                       gridPar1, dgridPar1,
                       &dataSize, dataField,
                       &out);
    }
  } else if (rank == 2) {
    if (use_fortran){
      FORTRAN_NAME(interpolate_2d_g)(&val_vec[0], &val_vec[1],
                                     gridDim,
                                     gridPar1, &dgridPar1,
                                     gridPar2, &dgridPar2,
                                     &dataSize, dataField,
                                     &out);
    } else {
      interpolate_2d_g(val_vec[0], val_vec[1],
                       gridDim,
                       gridPar1, dgridPar1,
                       gridPar2, dgridPar2,
                       &dataSize, dataField,
                       &out);
    }
  } else if (rank == 3) {
    if (use_fortran){
      FORTRAN_NAME(interpolate_3d_g)(&val_vec[0], &val_vec[1], &val_vec[2],
                                     gridDim,
                                     gridPar1, &dgridPar1,
                                     gridPar2, &dgridPar2,
                                     gridPar3, &dgridPar3,
                                     &dataSize, dataField,
                                     &out);
    } else {
      interpolate_3d_g(val_vec[0], val_vec[1], val_vec[2],
                       gridDim,
                       gridPar1, dgridPar1,
                       gridPar2, dgridPar2,
                       gridPar3, dgridPar3,
                       &dataSize, dataField,
                       &out);
    }
  } else if (rank == 4) {
    if (use_fortran){
      FORTRAN_NAME(interpolate_4d_g)(&val_vec[0], &val_vec[1], &val_vec[2],
                                     &val_vec[3],
                                     gridDim,
                                     gridPar1, &dgridPar1,
                                     gridPar2, &dgridPar2,
                                     gridPar3, &dgridPar3,
                                     gridPar4, &dgridPar4,
                                     &dataSize, dataField,
                                     &out);
    } else {
      interpolate_4d_g(val_vec[0], val_vec[1], val_vec[2],
                       val_vec[3],
                       gridDim,
                       gridPar1, dgridPar1,
                       gridPar2, dgridPar2,
                       gridPar3, dgridPar3,
                       gridPar4, dgridPar4,
                       &dataSize, dataField,
                       &out);
    }
  } else {
    if (use_fortran){
      FORTRAN_NAME(interpolate_5d_g)(&val_vec[0], &val_vec[1], &val_vec[2],
                                     &val_vec[3], &val_vec[4],
                                     gridDim,
                                     gridPar1, &dgridPar1,
                                     gridPar2, &dgridPar2,
                                     gridPar3, &dgridPar3,
                                     gridPar4, &dgridPar4,
                                     gridPar5, &dgridPar5,
                                     &dataSize, dataField,
                                     &out);
    } else {
      interpolate_5d_g(val_vec[0], val_vec[1], val_vec[2],
                       val_vec[3], val_vec[4],
                       gridDim,
                       gridPar1, dgridPar1,
                       gridPar2, dgridPar2,
                       gridPar3, dgridPar3,
                       gridPar4, dgridPar4,
                       gridPar5, dgridPar5,
                       &dataSize, dataField,
                       &out);
    }
  }
  return out;
}

std::vector<double> run_interp(const std::vector<std::vector<double>>& vals,
                               const InterpTable& table, bool use_fortran)
{
  std::vector<double> out(vals.size());
  for (std::size_t i = 0; i < vals.size(); i++) {
    out[i] = run_interp(vals[i], table, use_fortran);
  }
  return out;
}

double run_interp_3dz(const std::vector<double>& val_vec,
                      const InterpTable& table,
                      bool use_fortran)
{
  const int rank = table.ndim();
  if ((val_vec.size() != rank) && (rank != 3)){
    error("run_interp_3dz", "val_vec must have 3 parameters.");
  }

  const double* gridPar1 = table.gridPar_0Indexed(0);
  const double* gridPar2 = table.gridPar_0Indexed(1);
  const double* gridPar3 = table.gridPar_0Indexed(2);

  const double dgridPar1 = table.dgridPar_0Indexed(0);
  const double dgridPar3 = table.dgridPar_0Indexed(2);

  const gr_int64* gridDim = table.gridDim();
  const gr_int64 dataSize = table.dataSize();
  const double* dataField = table.gridField();

  const double zr = val_vec[1]; // z redshift

  // Calculate index for redshift dimension - intentionally kept 1-indexed
  const long long zindex = find_zindex(zr, rank, gridDim, gridPar2);
  const gr_int64 end_int = ((rank > 2) && (zindex == gridDim[1]));

  double out;
  if (use_fortran){
    FORTRAN_NAME(interpolate_3dz_g)(&val_vec[0], &zr, &val_vec[2],
                                    gridDim,
                                    gridPar1, &dgridPar1,
                                    gridPar2, &zindex,
                                    gridPar3, &dgridPar3,
                                    &dataSize, dataField,
                                    &end_int, &out);
  } else {
    interpolate_3dz_g(val_vec[0], zr, val_vec[2],
                      gridDim,
                      gridPar1, dgridPar1,
                      gridPar2, zindex,
                      gridPar3, dgridPar3,
                      &dataSize, dataField,
                      &end_int, &out);
  }
  return out;
}

std::vector<double> run_interp_3dz(const std::vector<std::vector<double>>& vals,
                                   const InterpTable& table, bool use_fortran)
{
  std::vector<double> out(vals.size());
  for (std::size_t i = 0; i < vals.size(); i++) {
    out[i] = run_interp_3dz(vals[i], table, use_fortran);
  }
  return out;
}

// uniform distribution on the interval (0, 1]
double uniform_dist_transform_(std::minstd_rand &prng){

  // sanity check to confirm that the largest value returned by generator is
  // representable (without any loss of precision) by a double (its <= 2^53)
  static_assert( (prng.max() <= 9007199254740992) & (prng.min() == 1),
                 "sanity check failed");

  return ( static_cast<double>(prng()) / static_cast<double>(prng.max()) );
}


struct ParamProps { double min_val; double max_val; int num_vals; };

InterpTable build_table(std::uint32_t seed,
                        const std::vector<ParamProps>& paramprop_vec)
{
  std::minstd_rand prng = std::minstd_rand(seed);
  std::vector<std::vector<double>> param_vals;
  std::size_t field_size = 1;

  for (const ParamProps& cur_param_prop : paramprop_vec) {
    const int dim_size = cur_param_prop.num_vals;

    // initialize cur_param_vals
    std::vector<double> cur_param_vals(dim_size);
    double dt = (cur_param_prop.max_val - cur_param_prop.min_val) / dim_size;
    for (int i = 0; i < dim_size; i++){
      cur_param_vals[i] = cur_param_prop.min_val + i * dt;
    }
    cur_param_vals[dim_size - 1] = cur_param_prop.max_val;
    param_vals.push_back(cur_param_vals);

    // update field_size
    field_size *= dim_size;
  }

  // now initialize field_data
  std::vector<double> field_data(field_size);
  for (int i = 0; i < field_size; i++) {
    field_data[i] = 100.0 * uniform_dist_transform_(prng) - 50.0;
  }

  return InterpTable(param_vals, field_data);
}

InterpTable get_3dz_table()
{
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

// this returns a vector of vectors appropriate interpolate function for cases where 1 or more
// inputs come from special_vals (the other inputs are taken from
// ordinary_vals)
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
std::vector<std::vector<double>> get_combinations_
(const std::vector<double>& special_vals,
 const std::vector<double>& ordinary_vals)
{
  std::size_t rank = special_vals.size();

  int choice[5] = {0, 0, 0, 0, 0};

  int out_count = pow(2, rank) - 1;
  std::vector<std::vector<double>> out;
  out.reserve(out_count - 1);

  while (out.size() < out_count) {
    std::vector<double> cur(rank);
    for (int i = 0; i < rank; i++){
      cur[i] = (choice[i] == 0) ? special_vals[i] : ordinary_vals[i];
    }
    out.push_back(cur);

    // now increment choice
    for (int i = 0; i < rank; i++){
      if (choice[i] == 0) { choice[i] = 1; break; }
      choice[i] = 0;
    }
  }

  return out;
}


void run_test(const InterpTable& table, const bool use_3dz = false)
{
  const int rank = table.ndim();

  auto compare_funcs = [&](const std::vector<std::vector<double>> &inputs)
    {
      std::vector<double> ref, actual;
      if (use_3dz) {
        ref = run_interp_3dz(inputs, table, true);
        actual = run_interp_3dz(inputs, table, false);
      } else {
        ref = run_interp(inputs, table, true);
        actual = run_interp(inputs, table, false);
      }

      if (false) {
        std::string ref_s = vec_to_string(ref);
        std::string actual_s = vec_to_string(actual);
        std::printf("   reference: %s\n", ref_s.c_str());
        std::printf("   actual: %s\n", actual_s.c_str());
      }
      compare_values(actual, ref, 0.0, 0.0, "");
    };

  // get useful values for test:
  std::vector<double> min_param_val(rank);
  std::vector<double> max_param_val(rank);
  std::vector<double> unaligned_inrange(rank);
  std::vector<double> unaligned_inrange_alt(rank);
  std::vector<double> half_min_param_val(rank);
  std::vector<double> double_max_param_val(rank);
  for (int i = 0; i < rank; i++){
    const double* cur_param_vals = table.gridPar_0Indexed(i);
    const int cur_param_len = table.gridDim()[i];

    min_param_val[i] = cur_param_vals[0];
    max_param_val[i] = cur_param_vals[cur_param_len-1];

    double delta = cur_param_vals[1] - cur_param_vals[0];

    unaligned_inrange[i] = cur_param_vals[0] + delta * 0.3;
    unaligned_inrange_alt[i] = cur_param_vals[0] + delta * 0.5;

    half_min_param_val[i] = min_param_val[i]/2;
    double_max_param_val[i] = max_param_val[i]*2;
  }

  std::printf(" -> comparing some in-range values\n");
  {
    std::vector<std::vector<double>> inputs
      = get_combinations_(unaligned_inrange, unaligned_inrange_alt);
    compare_funcs(inputs);
  }

  std::printf(" -> comparing cases with 1+ inputs on left grid edge\n");
  {
    std::vector<std::vector<double>> inputs
      = get_combinations_(min_param_val, unaligned_inrange);
    compare_funcs(inputs);
  }

  std::printf(" -> comparing cases with 1+ inputs on right grid edge\n");
  {
    std::vector<std::vector<double>> inputs
      = get_combinations_(max_param_val, unaligned_inrange);
    compare_funcs(inputs);
  }

  if (!use_3dz) {
    std::printf(" -> comparing cases with 1+ inputs less than min grid val\n");
    std::vector<std::vector<double>> inputs
      = get_combinations_(half_min_param_val, unaligned_inrange);
    compare_funcs(inputs);
  }

  std::printf(" -> comparing cases with 1+ inputs greater than min grid val\n");
  {
    std::vector<std::vector<double>> inputs
      = get_combinations_(double_max_param_val, unaligned_inrange);
    compare_funcs(inputs);
  }


  if (use_3dz) {
    std::printf(" -> some extra redshifts\n");

    gr_int64 n_gridded_z_vals = table.gridDim()[1];

    double temp[5] = {table.gridPar_0Indexed(1)[n_gridded_z_vals - 4],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 3],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 2],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 1],
                      table.gridPar_0Indexed(1)[n_gridded_z_vals - 1] * 1.01};

    std::vector<std::vector<double>> inputs;
    for (int i = 0; i < 5; i++) {
      inputs.push_back({unaligned_inrange[0], temp[i], unaligned_inrange[2]});
    }
    compare_funcs(inputs);
  }

}


int main(){

  const std::uint32_t seed = 342;


  std::printf("comparing interpolate_1d_g:\n");
  {
    InterpTable table = build_table(seed, {{-6.0, 6.0, 25}});
    run_test(table);
  }


  std::printf("\ncomparing interpolate_2d_g:\n");
  {
    InterpTable table = build_table(seed, {{-6.0, 6.0, 25}, {0.0, 10.0, 11}});
    run_test(table);
  }


  std::printf("\ncomparing interpolate_3d_g:\n");
  {
    InterpTable table = build_table(seed,
                                    {{-6.0, 6.0, 25}, {0.0, 10.0, 11},
                                     {-1.0, 0.0, 5}});
    run_test(table);
  }


  std::printf("\ncomparing interpolate_4d_g:\n");
  {
    InterpTable table = build_table(seed,
                                    {{-6.0, 6.0, 25}, {0.0, 10.0, 11},
                                     {-1.0, 0.0, 5}, {0.5, 5.5, 6}});
    run_test(table);
  }


  std::printf("\ncomparing interpolate_5d_g:\n");
  {
    InterpTable table = build_table(seed,
                                    {{-6.0, 6.0, 25}, {0.0, 10.0, 11},
                                     {-1.0, 0.0, 5}, {0.5, 5.5, 6},
                                     {-10.0, 0.0, 11}});
    run_test(table);
  }


  std::printf("\ncomparing interpolate_3dz_g:\n");
  {
    InterpTable table = get_3dz_table();
    run_test(table, true);
  }
  
  return 0;
}
