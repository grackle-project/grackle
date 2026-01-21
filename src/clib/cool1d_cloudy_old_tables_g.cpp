// See LICENSE file for license and copyright information

/// @file cool1d_cloudy_old_tables_g-cpp.C
/// @brief Declares signature of cool1d_cloudy_old_tables_g

// This file was initially generated automatically during conversion of the
// cool1d_cloudy_old_tables_g function from FORTRAN to C++

#include <cmath>
#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "fortran_func_wrappers.hpp"
#include "utils-cpp.hpp"

#include "cool1d_cloudy_old_tables_g.hpp"

void grackle::impl::cool1d_cloudy_old_tables_g(
    const double* rhoH, double* metallicity, const double* logtem, double* edot,
    double comp2, double dom, double zr, const gr_mask_type* itmask,
    chemistry_data* my_chemistry, cloudy_data cloudy_table, gr_float* density,
    gr_float* e_density, grackle_field_data* my_fields, IndexRange idx_range) {
  // General Arguments

  grackle::impl::View<gr_float***> d(density, idx_range.i_stop,
                                     my_fields->grid_dimension[1],
                                     my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> de(e_density, idx_range.i_stop,
                                      my_fields->grid_dimension[1],
                                      my_fields->grid_dimension[2]);

  // Locals

  int i;
  double inv_log10, log10_tCMB;
  double dclPar[GRACKLE_CLOUDY_TABLE_MAX_DIMENSION] = {};

  // Slice locals

  std::vector<double> log_Z(idx_range.i_stop);
  std::vector<double> e_frac(idx_range.i_stop);
  std::vector<double> log_e_frac(idx_range.i_stop);
  std::vector<double> cl_e_frac(idx_range.i_stop);
  std::vector<double> fh(idx_range.i_stop);
  std::vector<double> log_n_h(idx_range.i_stop);
  std::vector<double> log_cool(idx_range.i_stop);
  std::vector<double> log_cool_cmb(idx_range.i_stop);
  std::vector<double> log_heat(idx_range.i_stop);
  std::vector<double> edot_met(idx_range.i_stop);
  std::vector<double> log10tem(idx_range.i_stop);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  inv_log10 = 1. / std::log(10.);
  log10_tCMB = std::log10(comp2);

  // Calculate parameter value slopes

  dclPar[0] =
      (cloudy_table.grid_parameters[0][cloudy_table.grid_dimension[0] - 1] -
       cloudy_table.grid_parameters[0][0]) /
      (double)(cloudy_table.grid_dimension[0] - 1);
  if (cloudy_table.grid_rank > 1) {
    dclPar[1] =
        (cloudy_table.grid_parameters[1][cloudy_table.grid_dimension[1] - 1] -
         cloudy_table.grid_parameters[1][0]) /
        (double)(cloudy_table.grid_dimension[1] - 1);
  }
  if (cloudy_table.grid_rank > 2) {
    dclPar[2] =
        (cloudy_table.grid_parameters[2][cloudy_table.grid_dimension[2] - 1] -
         cloudy_table.grid_parameters[2][0]) /
        (double)(cloudy_table.grid_dimension[2] - 1);
  }
  if (cloudy_table.grid_rank > 3) {
    dclPar[3] =
        (cloudy_table.grid_parameters[3][cloudy_table.grid_dimension[3] - 1] -
         cloudy_table.grid_parameters[3][0]) /
        (double)(cloudy_table.grid_dimension[3] - 1);
  }
  if (cloudy_table.grid_rank > 4) {
    dclPar[4] =
        (cloudy_table.grid_parameters[4][cloudy_table.grid_dimension[4] - 1] -
         cloudy_table.grid_parameters[4][0]) /
        (double)(cloudy_table.grid_dimension[4] - 1);
  }

  for (i = idx_range.i_start; i <= idx_range.i_end; i++) {
    if (itmask[i] != MASK_FALSE) {
      log10tem[i] = logtem[i] * inv_log10;

      // Calcualte H mass fraction

      fh[i] = rhoH[i] / d(i, idx_range.j, idx_range.k);

      // Calculate proper log(n_H)

      if (cloudy_table.grid_rank > 1) {
        log_n_h[i] = std::log10(rhoH[i] * dom);
      }

      // Calculate metallicity

      if (cloudy_table.grid_rank > 2) {
        log_Z[i] = std::log10(metallicity[i]);
      }

      // Calculate electron fraction

      if (cloudy_table.grid_rank > 3) {
        e_frac[i] = 2. * de(i, idx_range.j, idx_range.k) /
                    (d(i, idx_range.j, idx_range.k) * (1. + fh[i]));
        // Make sure electron fraction is never above 1
        // which can give bad cooling/heating values when
        // extrapolating in the Cloudy data.
        log_e_frac[i] = std::min(std::log10(e_frac[i]), 0.0);

        // Get extra electrons contributed by metals

        cl_e_frac[i] =
            e_frac[i] *
            (1. + (2. * my_chemistry->cloudy_electron_fraction_factor *
                   metallicity[i] * fh[i]) /
                      (1. + fh[i]));
      }

      // Call interpolation functions to get heating/cooling

      // Interpolate over temperature.
      if (cloudy_table.grid_rank == 1) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
            log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1) &&
            ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_1d_g(
              log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density and temperature.
      } else if (cloudy_table.grid_rank == 2) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
            log_n_h[i], log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0],
            cloudy_table.grid_parameters[1], dclPar[1], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1) &&
            ((log10tem[i] - log10_tCMB) < 2.0f)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log_n_h[i], log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_2d_g(
              log_n_h[i], log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density, metallicity, and temperature.
      } else if (cloudy_table.grid_rank == 3) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_3d_g(
            log_n_h[i], log_Z[i], log10tem[i], cloudy_table.grid_dimension,
            cloudy_table.grid_parameters[0], dclPar[0],
            cloudy_table.grid_parameters[1], dclPar[1],
            cloudy_table.grid_parameters[2], dclPar[2], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1) &&
            ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_3d_g(
              log_n_h[i], log_Z[i], log10_tCMB, cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_3d_g(
              log_n_h[i], log_Z[i], log10tem[i], cloudy_table.grid_dimension,
              cloudy_table.grid_parameters[0], dclPar[0],
              cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density, metallicity, electron fraction, and
        // temperature.
      } else if (cloudy_table.grid_rank == 4) {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_4d_g(
            log_n_h[i], log_Z[i], log_e_frac[i], log10tem[i],
            cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
            dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
            cloudy_table.grid_parameters[2], dclPar[2],
            cloudy_table.grid_parameters[3], dclPar[3], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1) &&
            ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_4d_g(
              log_n_h[i], log_Z[i], log_e_frac[i], log10_tCMB,
              cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
              dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.grid_parameters[3], dclPar[3],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_4d_g(
              log_n_h[i], log_Z[i], log_e_frac[i], log10tem[i],
              cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
              dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.grid_parameters[3], dclPar[3],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }

        // Interpolate over density, metallicity, electron fraction, redshift,
        // and temperature.
      } else {
        log_cool[i] = grackle::impl::fortran_wrapper::interpolate_5d_g(
            log_n_h[i], log_Z[i], log_e_frac[i], zr, log10tem[i],
            cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
            dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
            cloudy_table.grid_parameters[2], dclPar[2],
            cloudy_table.grid_parameters[3], dclPar[3],
            cloudy_table.grid_parameters[4], dclPar[4], cloudy_table.data_size,
            cloudy_table.cooling_data);
        edot_met[i] = -std::pow(10., log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1) &&
            ((log10tem[i] - log10_tCMB) < 2.)) {
          log_cool_cmb[i] = grackle::impl::fortran_wrapper::interpolate_5d_g(
              log_n_h[i], log_Z[i], log_e_frac[i], zr, log10_tCMB,
              cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
              dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.grid_parameters[3], dclPar[3],
              cloudy_table.grid_parameters[4], dclPar[4],
              cloudy_table.data_size, cloudy_table.cooling_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1) {
          log_heat[i] = grackle::impl::fortran_wrapper::interpolate_5d_g(
              log_n_h[i], log_Z[i], log_e_frac[i], zr, log10tem[i],
              cloudy_table.grid_dimension, cloudy_table.grid_parameters[0],
              dclPar[0], cloudy_table.grid_parameters[1], dclPar[1],
              cloudy_table.grid_parameters[2], dclPar[2],
              cloudy_table.grid_parameters[3], dclPar[3],
              cloudy_table.grid_parameters[4], dclPar[4],
              cloudy_table.data_size, cloudy_table.heating_data);
          edot_met[i] = edot_met[i] + std::pow(10., log_heat[i]);
        }
      }

      if (cloudy_table.grid_rank > 3) {
        edot_met[i] = edot_met[i] * cl_e_frac[i];
      }

      edot[i] =
          edot[i] + (edot_met[i] * rhoH[i] * d(i, idx_range.j, idx_range.k));
    }
  }

  return;
}
