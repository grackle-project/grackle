#include <math.h>
#include <unistd.h>

#include "../utest_helpers.hpp"

extern "C" {
  #include "../../clib/interop/interop_funcs.h"

  // legacy version of the function
  extern void FORTRAN_NAME(calc_temp1d_cloudy_g)(
        const gr_float* d, const gr_float* metal, // 3D array
        const gr_float* e, // 3D array
        const double* rhoH, // 1D array
        const int* in, const int* jn, const int* kn,
        const int* is, const int* ie, const int* j, const int* k,
        double* tgas, double*mmw,
        const double* dom, const double* zr,
        const double* temstart, const double* temend, const double* gamma,
        const double* utem, const int* imetal,
        const long long* clGridRank,
        const long long* clGridDim,
        const double* clPar1,
        const double* clPar2,
        const double* clPar3,
        const long long* clDataSize,
        const double* clMMW,
        const int32_t* itmask);
}


struct calc_temp_outputs{
  std::vector<double> tgas;
  std::vector<double> mmw;
};

calc_temp_outputs run_test(DummyGrackleConfig& config,
                           bool use_fortran, bool slc_from_3D_arr = true){

  /* Calculate temperature units. */
  const double temperature_units = get_temperature_units(&(config.units));

  // initialize physical quantities!
  const gr_float Temp0 = 1.0e3;   // background temperature (in K)

  const std::size_t length = 60;
  std::vector<gr_float> density(length);
  std::vector<gr_float> metal_density(length);
  std::vector<gr_float> eint(length);
  std::vector<double> rhoH(length); // double is NOT a typo
  for (std::size_t i = 0; i < length; i++) {
    density[i] = gr_float(i+1);
    metal_density[i] = gr_float(config.chem_data.SolarMetalFractionByMass *
                                density[i]);
    rhoH[i] = config.chem_data.HydrogenFractionByMass * density[i];
    eint[i] = gr_float(Temp0 / temperature_units);
  }

  // prepare some other arguments
  //    -> when running in 3D, we aren't
  const int in = (!slc_from_3D_arr) ? int(length) : 5;
  const int jn = (!slc_from_3D_arr) ? 1           : 4;
  const int kn = (!slc_from_3D_arr) ? 1           : 3;

  if (std::size_t(in * jn * kn) > length) { error("something is wrong\n"); }

  const int j = jn; // remember, these are currently 1-indexed
  const int k = kn;

  const int is = 0; // not a typo! (even though other indices are 1-indexed)
  const int ie = in - 1; // not a typo!

  const double aye = config.units.a_value; // expansion factor (in code units)
  const double dom = config.units.density_units * pow(aye,3) / mh;
  const double zr = 1.0/(aye * config.units.a_units) - 1.;
  const int imetal = 1;

  // prepare the output arguments
  std::vector<double> tgas(in);
  std::vector<double> mmw(in);

  const std::vector<int32_t> itmask(std::size_t(in), int32_t(1));

  if (use_fortran) {
    FORTRAN_NAME(calc_temp1d_cloudy_g)
      (density.data(), metal_density.data(), eint.data(),
       rhoH.data(),
       &in, &jn, &kn,  &is, &ie, &j, &k,
       tgas.data(), mmw.data(),
       &dom, &zr,
       &config.chem_data.TemperatureStart, &config.chem_data.TemperatureEnd,
       &config.chem_data.Gamma, &temperature_units, &imetal,
       &config.chem_rates.cloudy_primordial.grid_rank, // clGridRank
       config.chem_rates.cloudy_primordial.grid_dimension, // clGridDim
       config.chem_rates.cloudy_primordial.grid_parameters[0], // clPar1
       config.chem_rates.cloudy_primordial.grid_parameters[1], // clPar2
       config.chem_rates.cloudy_primordial.grid_parameters[2], // clPar3
       &config.chem_rates.cloudy_primordial.data_size, // clDataSize
       config.chem_rates.cloudy_primordial.mmw_data, // clMMW
       itmask.data());

  } else {
    calc_temp1d_cloudy_g
      (density.data(), metal_density.data(), eint.data(),
       rhoH.data(),
       in, jn, kn,  is, ie, j, k,
       tgas.data(), mmw.data(),
       dom, zr,
       config.chem_data.TemperatureStart, config.chem_data.TemperatureEnd,
       config.chem_data.Gamma, temperature_units, imetal,
       config.chem_rates.cloudy_primordial.grid_rank, // clGridRank
       config.chem_rates.cloudy_primordial.grid_dimension, // clGridDim
       config.chem_rates.cloudy_primordial.grid_parameters[0], // clPar1
       config.chem_rates.cloudy_primordial.grid_parameters[1], // clPar2
       config.chem_rates.cloudy_primordial.grid_parameters[2], // clPar3
       config.chem_rates.cloudy_primordial.data_size, // clDataSize
       config.chem_rates.cloudy_primordial.mmw_data, // clMMW
       itmask.data());
  }

  return {tgas, mmw};
}

int main(void){
  for (int n_tab_dims = 1; n_tab_dims < 4; n_tab_dims++) {
    printf("\nConsidering a %dD table of MMW values\n", n_tab_dims);

    std::vector<double> z_vals {0.0};
    if (n_tab_dims == 3) {
      z_vals.push_back(0.13242);
      z_vals.push_back(0.5);
      z_vals.push_back(10);

      DummyGrackleConfig conf(n_tab_dims,0.0);
      double* clPar2 = conf.chem_rates.cloudy_primordial.grid_parameters[1];
      long long* clGridDim = conf.chem_rates.cloudy_primordial.grid_dimension;
      z_vals.push_back(clPar2[clGridDim[1]-3]);
      z_vals.push_back(clPar2[clGridDim[1]-2]);
      z_vals.push_back(clPar2[clGridDim[1]-1]);
      z_vals.push_back(clPar2[clGridDim[1]-1] * 1.01);
    }

    for (const auto& z_val : z_vals) {
      DummyGrackleConfig config(n_tab_dims,z_val);

      for (bool slc_from_3D_arr : {false, true}) {
        const char* descr = (slc_from_3D_arr) ? "slice from 3D arr" : "1D arr";
        printf("-> z = %g, comparing %s\n", z_val, descr);

        //printf("using the c version:\n"); fflush(stdout);
        calc_temp_outputs actual = run_test(config, false, slc_from_3D_arr);

        //printf("using the fortran version:\n"); fflush(stdout);
        calc_temp_outputs reference = run_test(config, true, slc_from_3D_arr);

        compare_values(actual.mmw, reference.mmw, 0.0, 0.0,
                       "**Error during comparison of mmw**");
        compare_values(actual.tgas, reference.tgas, 0.0, 0.0,
                       "**Error during comparison of tgas**");
      }
    }
  }

  return 0;
}
