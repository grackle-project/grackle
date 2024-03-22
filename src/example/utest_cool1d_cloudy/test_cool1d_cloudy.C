#include <math.h>
#include <unistd.h>

#include "../utest_helpers.hpp"

extern "C" {
  #include "../../clib/interop/interop_funcs.h"

  extern void FORTRAN_NAME(cool1d_cloudy_g)(
        const gr_float* d, // 3D arrays
        const double* rhoH, const gr_float* metallicity, // 1D array
        const int* in, const int* jn, const int* kn,
        const int* is, const int* ie, const int* j, const int* k,
        const double* logtem, double* edot, // 1D array
        const double* comp2, const double* dom, const double* zr,
        const int* icmbTfloor, const int* iClHeat, const int* iZscale,
        const long long* clGridRank,
        const long long* clGridDim,
        const double* clPar1, const double* clPar2, const double* clPar3,
        const long long* clDataSize,
        const double* clCooling, const double* clHeating,
        const int32_t* itmask);
}

std::vector<double> run_test(DummyGrackleConfig& config,
                             bool use_metal_table,
                             bool use_fortran,
                             bool slc_from_3D_arr = true){

  /* Calculate temperature units. */
  const double temperature_units = get_temperature_units(&(config.units));

  // initialize physical quantities!
  const gr_float Temp0 = 1.0e3;   // background temperature (in K)

  const double z_solar = config.chem_data.SolarMetalFractionByMass;

  const std::size_t length = 60;
  std::vector<gr_float> density(length);
  std::vector<gr_float> metal_density(length);
  std::vector<gr_float> eint(length);
  std::vector<double> rhoH(length); // double is NOT a typo
  for (std::size_t i = 0; i < length; i++) {
    density[i] = gr_float(i+1);
    metal_density[i] = gr_float(z_solar * density[i]);
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

  // calculate temperature (as part of the setup)
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

  // prepare a couple other arguments!
  std::vector<double> metallicity(length); // double is NOT a typo
  for (std::size_t i = 0; i < length; i++) {
    metallicity[i] = metal_density[i] / density[i] / z_solar;
  }

  // these use natural logarithms (intentional!)
  const double logtem0 = log(config.chem_data.TemperatureStart);
  const double logtem9 = log(config.chem_data.TemperatureEnd);

  std::vector<double> logtem(length); // double is NOT a typo
  for (std::size_t i = 0; i < length; i++) {
    if (tgas[i] <= 0) { logtem[i] = logtem0; } // tgas wasn't computed here
    logtem[i] = fmin(fmax(log(tgas[i]), logtem0), logtem9);
  }

  const int iClHeat = config.chem_data.UVbackground;

  const cloudy_data& cloudy_data = (use_metal_table) ?
    config.chem_rates.cloudy_metal : config.chem_rates.cloudy_primordial;

  const int iZscale = (use_metal_table) ? 1 : 0;
  const int icmbTfloor = (use_metal_table) ?
    config.chem_data.cmb_temperature_floor : 0;

  // things to vary (in the future):
  // metal cooling: icmbTfloor (when using metal cooling)
  // in both cases, vary iClHeat

  // compton cooling coefficient
  const double comp2 = 2.73 * (1.0 + zr);

  std::vector<double> edot(length, 0.0); // all outputs are set to 0
  // now we actually compute the cooling time
  if (use_fortran) {
    FORTRAN_NAME(cool1d_cloudy_g)(
        density.data(), // 3D arrays
        rhoH.data(), metallicity.data(), // 1D array
        &in, &jn, &kn, &is, &ie, &j, &k,
        logtem.data(), edot.data(), // 1D array
        &comp2, &dom, &zr,
        &icmbTfloor, &iClHeat, &iZscale,
        &cloudy_data.grid_rank, // clGridRank
        cloudy_data.grid_dimension, // clGridDim
        cloudy_data.grid_parameters[0], // clPar1
        cloudy_data.grid_parameters[1], // clPar2
        cloudy_data.grid_parameters[2], // clPar3
        &cloudy_data.data_size, // clDataSize
        cloudy_data.cooling_data, // clCooling
        cloudy_data.heating_data, // clHeating
        itmask.data());
  } else {
    cool1d_cloudy_g(
        density.data(), // 3D arrays
        rhoH.data(), metallicity.data(), // 1D array
        in, jn, kn, is, ie, j, k,
        logtem.data(), edot.data(), // 1D array
        comp2, dom, zr,
        icmbTfloor, iClHeat, iZscale,
        cloudy_data.grid_rank, // clGridRank
        cloudy_data.grid_dimension, // clGridDim
        cloudy_data.grid_parameters[0], // clPar1
        cloudy_data.grid_parameters[1], // clPar2
        cloudy_data.grid_parameters[2], // clPar3
        cloudy_data.data_size, // clDataSize
        cloudy_data.cooling_data, // clCooling
        cloudy_data.heating_data, // clHeating
        itmask.data());
  }

  return edot;
}


int main(void){
  for (int n_tab_dims = 1; n_tab_dims < 4; n_tab_dims++) {
    printf("\nConsidering a %dD table of values\n", n_tab_dims);

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
      for (int j = 0; j < 3; j++){
        std::string descr;
        bool use_metal_table;
        bool cmb_temperature_floor;

        // there's no reason to apply the cmb_temperature_floor when
        // primordial chemistry is used

        if (j == 0) {
          use_metal_table = true;
          descr = "metals, cmb floor";
          cmb_temperature_floor = true;
        } else if (j == 1) {
          use_metal_table = true;
          descr = "metals, no cmb floor";
          cmb_temperature_floor = false;
        } else if (j == 2) {
          use_metal_table = false;
          descr = "primordial, no cmb floor";
          cmb_temperature_floor = true;
        } else {
          error("SOMETHING IS WRONG");
        }

        for (bool use_UVbackground : {false, true}) {
          std::string UV_descr;
          if (!use_UVbackground) {
            UV_descr = "cool-only";
          } else if (n_tab_dims == 3) {
            UV_descr = "heat & cool";
          } else {
            continue; // can't model UV background 
          }

          DummyGrackleConfig config(n_tab_dims,z_val, use_UVbackground,
                                    cmb_temperature_floor);

          for (bool slc_from_3D_arr : {false, true}) {
            const char* descr2 = (slc_from_3D_arr) ?
              "slice from 3D arr" : "1D arr";
            printf("-> z = %g, %s, %s, comparing %s\n",
                   z_val, descr.c_str(), UV_descr.c_str(), descr2);

            std::vector<double> actual = run_test(config, use_metal_table,
                                                  false, slc_from_3D_arr);
            std::vector<double> reference = run_test(config, use_metal_table,
                                                     true, slc_from_3D_arr);

            //printf("%.15e, %.15e\n", actual[0], reference[0]);

            compare_values(actual, reference, 0.0, 0.0,
                           "**Error during comparison of edot**");
          }
        }
      }
    }
  }
}
