#include <stdlib.h> // abort
#include <stdio.h>  // stderr, vfprintf
#include <stdarg.h> // va_list, va_start, va_end

#include <string>
#include <vector>

#define mh     1.67262171e-24

#ifndef OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
#endif

extern "C" {
  #include <grackle.h>
  #include "../clib/grackle_macros.h"
}

#undef max
#undef min

[[noreturn]] void error(const char *fmt, ...){

  // parse argument list
  va_list arg_list;
  va_start(arg_list, fmt);         // access variable arguments after fmt
  vfprintf(stderr, fmt, arg_list); // print to stderr
  va_end(arg_list);                // cleanup variable arguments
  abort();                         // exit program with nonzero exit code
}

inline std::string vec_to_string(const std::vector<double>& vec) {
  std::string out = "{";

  std::size_t len = vec.size();

  std::size_t pause_start;
  std::size_t pause_stop;

  if (len > 30){
    pause_start = 3;
    pause_stop = len - 3;
  } else {
    pause_start = len *2;
    pause_stop = pause_start;
  }

  for (std::size_t i = 0; i < len; i++) {
    if ((i > pause_start) && (i < pause_stop)) { continue; }

    if (i == pause_stop) {
      out += ", ... ";
    } else if (i != 0) {
      out += ", ";
    }

    char buf[30];
    sprintf(buf, "%g", vec[i]);
    out += buf;
  }
  return out + "}";
}

inline void compare_values(const std::vector<double>& actual,
                           const std::vector<double>& desired,
                           double rtol = 0.0, double atol = 0.0,
                           std::string err_msg = "")
{
  if (actual.size() != desired.size()){
    error("the compared arrays have different lengths\n");
  }

  std::size_t num_mismatches = 0;
  double max_absDiff = 0.0;
  std::size_t max_absDiff_ind = 0;
  double max_relDiff = 0.0;
  std::size_t max_relDiff_ind = 0;
  bool has_nan_mismatch = false;

  for (std::size_t i = 0; i < actual.size(); i++) {
    double cur_absDiff = fabs(actual[i]-desired[i]);

    bool isnan_actual = isnan(actual[i]);
    bool isnan_desired = isnan(desired[i]);

    if ( (cur_absDiff > (atol + rtol * fabs(desired[i]))) ||
         (isnan_actual != isnan_desired) ) {
      
      num_mismatches++;
      if (isnan_actual != isnan_desired){
        has_nan_mismatch = true;
        max_absDiff = NAN;
        max_absDiff_ind = i;
        max_relDiff = NAN;
        max_relDiff_ind = i;
      } else if (!has_nan_mismatch) {
        if (cur_absDiff > max_absDiff){
          max_absDiff = cur_absDiff;
          max_absDiff_ind = i;
        }

        if ( cur_absDiff > (max_relDiff * fabs(desired[i])) ) {
          max_relDiff = cur_absDiff / fabs(desired[i]);
          max_relDiff_ind = i;
        }
      }
    }
  }

  if (num_mismatches == 0) { return; }

  std::string actual_vec_str = vec_to_string(actual);
  std::string ref_vec_str = vec_to_string(desired);

  error
    (("arrays are unequal for the tolerance: rtol = %g, atol = %g\n"
      "%s\n" // custom error message
      "Mismatched elements: %d / %d\n"
      "Max absolute difference: %g,   ind = %d, actual = %g, reference = %g\n"
      "Max relative difference: %g,   ind = %d, actual = %g, reference = %g\n"
      "actual: %s\n"
      "desired: %s\n"),
     rtol, atol,
     err_msg.c_str(),
     (int)num_mismatches, (int)actual.size(),
     max_absDiff, (int)max_absDiff_ind, actual[max_absDiff_ind],
     desired[max_absDiff_ind],
     max_relDiff, (int)max_relDiff_ind, actual[max_relDiff_ind],
     desired[max_relDiff_ind],
     actual_vec_str.c_str(), ref_vec_str.c_str());
}

inline void cut_down_to_1D_table(cloudy_data* ptr){
  if (ptr->grid_rank != 2) { error("can currently only cut down 2D to 1D"); }

  // table over rho and Temperature
  // -> I confirmed from the mmw table that temperature axis is contiguous
  long long num_rho_vals = ptr->grid_dimension[0];
  long long num_T_vals = ptr->grid_dimension[1];

  cloudy_data newObj;
  newObj.grid_rank = 1;

  // Dimension of dataset
  newObj.grid_dimension
    = (long long *)malloc(sizeof(long long) * newObj.grid_rank);
  newObj.grid_dimension[0] = num_T_vals;
  GRACKLE_FREE(ptr->grid_dimension);

  // Dataset parameter values.
  newObj.grid_parameters
    = (double**)malloc(sizeof(double**) * 1);
  newObj.grid_parameters[0] = ptr->grid_parameters[1];
  GRACKLE_FREE(ptr->grid_parameters[0]);
  GRACKLE_FREE(ptr->grid_parameters);

  // since Temperature is the fast-index, we are going to just reuse the
  // original pointers (even though they will have a bunch of unused data)
  newObj.heating_data = ptr->heating_data;
  newObj.cooling_data = ptr->cooling_data;
  newObj.mmw_data = ptr->mmw_data;
  newObj.data_size = num_T_vals;

  (*ptr) = newObj;
}

struct DummyGrackleConfig{
  // the central purpose here is to hold grackle configuration options for
  // tabulated solver
  //
  // do NOT try to use this as a general purpose C++ interface

  chemistry_data chem_data;
  chemistry_data_storage chem_rates;
  code_units units;

  DummyGrackleConfig(int n_tab_dims, double radiation_redshit,
                     bool UVbackground = false,
                     bool cmb_temperature_floor = true) {
    // radiation_redshift is meaningless when n_tab_dims isn't 3
    
    // setup units!
    code_units my_units;
    my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
    my_units.density_units = mh;                  // = mass_H/cm^3 (in g/cm^3)
    my_units.length_units = 1.0 / 3.24077929e-25; // = 1 Mpc (in cm)
    my_units.time_units = 3.15576e13;             // = 1 Myr (in seconds)
    my_units.a_units = 1.0; // units for the expansion factor
    // Set expansion factor to 1 for non-cosmological simulation.
    my_units.a_value = 1. / (1. + radiation_redshit) / my_units.a_units;
    set_velocity_units(&my_units);

    // setup runtime parameters
    chemistry_data my_chem;
    local_initialize_chemistry_parameters(&my_chem);
    my_chem.use_grackle = 1;            // chemistry on
    my_chem.with_radiative_cooling = 1; // cooling on
    my_chem.primordial_chemistry = 0;
    my_chem.metal_cooling = 1;         // metal cooling on
    my_chem.UVbackground = (UVbackground) ? 1 : 0;          // UV background on
    my_chem.cmb_temperature_floor = (cmb_temperature_floor) ? 1 : 0;

    if (UVbackground && (n_tab_dims != 3)) {
      error("can't enable UVbackground without the redshift dependence");
    }

    if ((n_tab_dims <= 0) || (n_tab_dims > 3)) {
      error("n_tab_dims must lie be 1, 2, or 3\n");
    } else if (n_tab_dims <= 2) {
      my_chem.grackle_data_file
        = "../../grackle_data_files/input/CloudyData_noUVB.h5";
    } else {
      my_chem.grackle_data_file =
        "../../grackle_data_files/input/CloudyData_UVB=HM2012.h5";
    }

    chemistry_data_storage my_rates;
    local_initialize_chemistry_data(&my_chem, &my_rates, &my_units);

    if (n_tab_dims == 1) {
      cut_down_to_1D_table(&(my_rates.cloudy_metal));
      cut_down_to_1D_table(&(my_rates.cloudy_primordial));
      //error("Can't currently support n_tab_dims == 1\n");
    }

    this->chem_data = my_chem;
    this->chem_rates = my_rates;
    this->units = my_units;
  }

  ~DummyGrackleConfig() {
    local_free_chemistry_data(&(this->chem_data), &(this->chem_rates));
    // don't deallocate (this->chem_data).grackle_data_file, that member is
    // assigned a string-literal (it's not heap-allocated)
  }

  DummyGrackleConfig(const DummyGrackleConfig&) = delete;
  DummyGrackleConfig(DummyGrackleConfig&&) = delete;
  DummyGrackleConfig& operator=(const DummyGrackleConfig&) = delete;
  DummyGrackleConfig& operator=(DummyGrackleConfig&&) = delete;

};
