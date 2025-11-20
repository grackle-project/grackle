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
#include "utils-cpp.hpp"

#include "cool1d_cloudy_old_tables_g.hpp"

void grackle::impl::cool1d_cloudy_old_tables_g(
  double* rhoH, double* metallicity, double* logtem, double* edot,
  double comp2, double dom, double zr, gr_mask_type* itmask,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range
)
{

  // General Arguments

  grackle::impl::View<gr_float***> d(my_fields->density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);
  grackle::impl::View<gr_float***> de(my_fields->e_density, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // Locals

  int i, q;
  double inv_log10, log10_tCMB;
  std::vector<double> dclPar(my_rates->cloudy_metal.grid_rank);

  // Slice locals

  std::vector<double> log_Z(my_fields->grid_dimension[0]);
  std::vector<double> e_frac(my_fields->grid_dimension[0]);
  std::vector<double> log_e_frac(my_fields->grid_dimension[0]);
  std::vector<double> cl_e_frac(my_fields->grid_dimension[0]);
  std::vector<double> fh(my_fields->grid_dimension[0]);
  std::vector<double> log_n_h(my_fields->grid_dimension[0]);
  std::vector<double> log_cool(my_fields->grid_dimension[0]);
  std::vector<double> log_cool_cmb(my_fields->grid_dimension[0]);
  std::vector<double> log_heat(my_fields->grid_dimension[0]);
  std::vector<double> edot_met(my_fields->grid_dimension[0]);
  std::vector<double> log10tem(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================

  inv_log10 = 1. / std::log(10.);
  log10_tCMB = std::log10(comp2);

  // Calculate parameter value slopes

  dclPar[0] = (my_rates->cloudy_metal.grid_parameters[0][my_rates->cloudy_metal.grid_dimension[0]-1] - my_rates->cloudy_metal.grid_parameters[0][0]) /
       (double)(my_rates->cloudy_metal.grid_dimension[0] - 1 );
  if (my_rates->cloudy_metal.grid_rank > 1)  {
    dclPar[1] = (my_rates->cloudy_metal.grid_parameters[1][my_rates->cloudy_metal.grid_dimension[1]-1] - my_rates->cloudy_metal.grid_parameters[1][0]) /
         (double)(my_rates->cloudy_metal.grid_dimension[1] - 1 );
  }
  if (my_rates->cloudy_metal.grid_rank > 2)  {
    dclPar[2] = (my_rates->cloudy_metal.grid_parameters[2][my_rates->cloudy_metal.grid_dimension[2]-1] - my_rates->cloudy_metal.grid_parameters[2][0]) /
         (double)(my_rates->cloudy_metal.grid_dimension[2] - 1 );
  }
  if (my_rates->cloudy_metal.grid_rank > 3)  {
    dclPar[3] = (my_rates->cloudy_metal.grid_parameters[3][my_rates->cloudy_metal.grid_dimension[3]-1] - my_rates->cloudy_metal.grid_parameters[3][0]) /
         (double)(my_rates->cloudy_metal.grid_dimension[3] - 1 );
  }
  if (my_rates->cloudy_metal.grid_rank > 4)  {
    dclPar[4] = (my_rates->cloudy_metal.grid_parameters[4][my_rates->cloudy_metal.grid_dimension[4]-1] - my_rates->cloudy_metal.grid_parameters[4][0]) /
         (double)(my_rates->cloudy_metal.grid_dimension[4] - 1 );
  }

  for (i = idx_range.i_start; i<=idx_range.i_end ; i++) {
    if ( itmask[i] != MASK_FALSE )  {

      log10tem[i] = logtem[i] * inv_log10;

      // Calcualte H mass fraction

      fh[i] = rhoH[i] / d(i,idx_range.jp1-1,idx_range.kp1-1);

      // Calculate proper log(n_H)

      if (my_rates->cloudy_metal.grid_rank > 1)  {

        log_n_h[i] = std::log10(rhoH[i] * dom);

      }

      // Calculate metallicity

      if (my_rates->cloudy_metal.grid_rank > 2)  {

        log_Z[i] = std::log10(metallicity[i]);

      }

      // Calculate electron fraction
            
      if (my_rates->cloudy_metal.grid_rank > 3)  {

        e_frac[i] = 2. * de(i,idx_range.jp1-1,idx_range.kp1-1) /
             (d(i,idx_range.jp1-1,idx_range.kp1-1) * (1. + fh[i]));
        // Make sure electron fraction is never above 1
        // which can give bad cooling/heating values when
        // extrapolating in the Cloudy data.
        log_e_frac[i] = std::min(std::log10(e_frac[i]), 0.0);

        // Get extra electrons contributed by metals

        cl_e_frac[i] = e_frac[i] *
             (1. + (2. * my_chemistry->cloudy_electron_fraction_factor * metallicity[i] *
             fh[i]) / (1. + fh[i]));

      }

      // Call interpolation functions to get heating/cooling

      // Interpolate over temperature.
      if (my_rates->cloudy_metal.grid_rank == 1)  {
         FORTRAN_NAME(interpolate_1d_g)(
             &log10tem[i], my_rates->cloudy_metal.grid_dimension, my_rates->cloudy_metal.grid_parameters[0],
             &dclPar[0], &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool[i]);
        edot_met[i] = -std::pow(10.,log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1)  && 
             ((log10tem[i] - log10_tCMB) < 2.))  {
           FORTRAN_NAME(interpolate_1d_g)(
               &log10_tCMB, my_rates->cloudy_metal.grid_dimension, my_rates->cloudy_metal.grid_parameters[0],
               &dclPar[0], &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data,
               &log_cool_cmb[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1)  {
           FORTRAN_NAME(interpolate_1d_g)(
               &log10tem[i], my_rates->cloudy_metal.grid_dimension, my_rates->cloudy_metal.grid_parameters[0],
               &dclPar[0], &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.heating_data,
               &log_heat[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_heat[i]);
        }

        // Interpolate over density and temperature.
      } else if (my_rates->cloudy_metal.grid_rank == 2)  {
         FORTRAN_NAME(interpolate_2d_g)(
             &log_n_h[i], &log10tem[i], my_rates->cloudy_metal.grid_dimension,
             my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
             &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool[i]);
        edot_met[i] = -std::pow(10.,log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1)  && 
             ((log10tem[i] - log10_tCMB) < 2.0f))  {
           FORTRAN_NAME(interpolate_2d_g)(
               &log_n_h[i], &log10_tCMB, my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool_cmb[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1)  {
           FORTRAN_NAME(interpolate_2d_g)(
                  &log_n_h[i], &log10tem[i], my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.heating_data, &log_heat[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_heat[i]);
        }

        // Interpolate over density, metallicity, and temperature.
      } else if (my_rates->cloudy_metal.grid_rank == 3)  {
         FORTRAN_NAME(interpolate_3d_g)(
             &log_n_h[i], &log_Z[i], &log10tem[i],
             my_rates->cloudy_metal.grid_dimension,
             my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
             my_rates->cloudy_metal.grid_parameters[2], &dclPar[2],
             &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool[i]);
        edot_met[i] = -std::pow(10.,log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1)  && 
             ((log10tem[i] - log10_tCMB) < 2.))  {
           FORTRAN_NAME(interpolate_3d_g)(
               &log_n_h[i], &log_Z[i], &log10_tCMB,
               my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool_cmb[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1)  {
           FORTRAN_NAME(interpolate_3d_g)(
               &log_n_h[i], &log_Z[i], &log10tem[i],
               my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.heating_data, &log_heat[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_heat[i]);
        }

        // Interpolate over density, metallicity, electron fraction, and temperature.
      } else if (my_rates->cloudy_metal.grid_rank == 4)  {
         FORTRAN_NAME(interpolate_4d_g)(
             &log_n_h[i], &log_Z[i],
             &log_e_frac[i], &log10tem[i], my_rates->cloudy_metal.grid_dimension,
             my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
             my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
             &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool[i]);
        edot_met[i] = -std::pow(10.,log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1)  && 
             ((log10tem[i] - log10_tCMB) < 2.))  {
           FORTRAN_NAME(interpolate_4d_g)(
               &log_n_h[i], &log_Z[i],
               &log_e_frac[i], &log10_tCMB, my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool_cmb[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1)  {
           FORTRAN_NAME(interpolate_4d_g)(
               &log_n_h[i], &log_Z[i],
               &log_e_frac[i], &log10tem[i], my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.heating_data, &log_heat[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_heat[i]);
        }

        // Interpolate over density, metallicity, electron fraction, redshift,
        // and temperature.
      } else {
         FORTRAN_NAME(interpolate_5d_g)(
             &log_n_h[i], &log_Z[i],
             &log_e_frac[i], &zr, &log10tem[i], my_rates->cloudy_metal.grid_dimension,
             my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
             my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
             my_rates->cloudy_metal.grid_parameters[4], &dclPar[4],
             &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool[i]);
        edot_met[i] = -std::pow(10.,log_cool[i]);

        // Ignore CMB term if T >> T_CMB
        if ((my_chemistry->cmb_temperature_floor == 1)  && 
             ((log10tem[i] - log10_tCMB) < 2.))  {
           FORTRAN_NAME(interpolate_5d_g)(
               &log_n_h[i], &log_Z[i],
               &log_e_frac[i], &zr, &log10_tCMB, my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
               my_rates->cloudy_metal.grid_parameters[4], &dclPar[4],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, &log_cool_cmb[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_cool_cmb[i]);
        }

        if (my_chemistry->UVbackground == 1)  {
           FORTRAN_NAME(interpolate_5d_g)(
               &log_n_h[i], &log_Z[i],
               &log_e_frac[i], &zr, &log10tem[i], my_rates->cloudy_metal.grid_dimension,
               my_rates->cloudy_metal.grid_parameters[0], &dclPar[0], my_rates->cloudy_metal.grid_parameters[1], &dclPar[1],
               my_rates->cloudy_metal.grid_parameters[2], &dclPar[2], my_rates->cloudy_metal.grid_parameters[3], &dclPar[3],
               my_rates->cloudy_metal.grid_parameters[4], &dclPar[4],
               &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.heating_data, &log_heat[i]);
          edot_met[i] = edot_met[i] + std::pow(10.,log_heat[i]);
        }

      }

      if (my_rates->cloudy_metal.grid_rank > 3)  {
        edot_met[i] = edot_met[i] * cl_e_frac[i];
      }

      edot[i] = edot[i] + (edot_met[i] * rhoH[i] * d(i,idx_range.jp1-1,idx_range.kp1-1));

    }
  }

  return;
}
