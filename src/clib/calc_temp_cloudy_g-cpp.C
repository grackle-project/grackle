// See LICENSE file for license and copyright information

/// @file calc_temp_cloudy_g-cpp.C
/// @brief Declares signature of calc_temp_cloudy_g

// This file was initially generated automatically during conversion of the
// calc_temp_cloudy_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "utils-cpp.hpp"

#include "calc_temp_cloudy_g-cpp.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void calc_temp_cloudy_g(
  gr_float* temperature_data_, int* imetal, chemistry_data* my_chemistry,
  chemistry_data_storage* my_rates, grackle_field_data* my_fields,
  InternalGrUnits internalu
)
{

  grackle::impl::View<gr_float***> d(
    my_fields->density,
    my_fields->grid_dimension[0],
    my_fields->grid_dimension[1],
    my_fields->grid_dimension[2]
  );
  grackle::impl::View<gr_float***> temperature(
    temperature_data_,
    my_fields->grid_dimension[0],
    my_fields->grid_dimension[1],
    my_fields->grid_dimension[2]
  );

  const double mh_local_var = mh_grflt;

  // Locals

  int i, j, k;
  int t, dj, dk;
  double dom, zr;
  double dbase1, tbase1, xbase1;
  double factor;

  // row temporaries

  std::vector<double> tgas(my_fields->grid_dimension[0]);
  std::vector<double> rhoH(my_fields->grid_dimension[0]);
  std::vector<double> mmw(my_fields->grid_dimension[0]);

  // Iteration mask

  std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
  
  // Set units

  dom      = internalu.urho*(std::pow(internalu.a_value,3))/mh_local_var;
  tbase1   = internalu.tbase1;
  xbase1   = internalu.uxyz/(internalu.a_value*internalu.a_units);    // uxyz is [x]*a      = [x]*[a]*a'        '
  dbase1   = internalu.urho*std::pow((internalu.a_value*internalu.a_units),3); // urho is [dens]/a^3 = [dens]/([a]*a')^3 '
  zr       = 1./(internalu.a_value*internalu.a_units) - 1.;

  // Convert densities from comoving to proper

  if (internalu.extfields_in_comoving == 1)  {
    factor = std::pow(internalu.a_value,(-3));
     FORTRAN_NAME(scale_fields_table_g)(d.data(), my_fields->metal_density,
                   &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                   &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], imetal, &factor);

  }

  // Loop over zones, and do an entire i-column in one go

  dk = my_fields->grid_end[2] - my_fields->grid_start[2] + 1;
  dj = my_fields->grid_end[1] - my_fields->grid_start[1] + 1;

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: !$omp parallel do schedule(runtime) private(
  //_// PORT: !$omp&  i, j, k, 
  //_// PORT: !$omp&  tgas, rhoH, mmw,
  //_// PORT: !$omp&  itmask )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    //_// TODO: move relevant variable declarations to here to replace OMP private
    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (t = 0; t<=(dk * dj - 1); t++) {
      k = t/dj      + my_fields->grid_start[2]+1;
      j = grackle::impl::mod(t,dj) + my_fields->grid_start[1]+1;

      // Initialize iteration mask to true for all cells.

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        itmask[i-1] = MASK_TRUE;

        rhoH[i-1] = my_chemistry->HydrogenFractionByMass * d(i-1,j-1,k-1);
      }

      // Calculate temperature and mean molecular weight

       FORTRAN_NAME(calc_temp1d_cloudy_g)(d.data(), my_fields->metal_density, my_fields->internal_energy, rhoH.data(),
           &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &j, &k,
           tgas.data(), mmw.data(), &dom, &zr,
           &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
           &my_chemistry->Gamma, &internalu.utem, imetal,
           &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
           my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2],
           &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.mmw_data,
           itmask.data());

      // Copy slice values into field array

      for (i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        temperature(i-1,j-1,k-1) = tgas[i-1];
      }

    }
  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from proper

  if (internalu.extfields_in_comoving == 1)  {
    factor = std::pow(internalu.a_value,3);
     FORTRAN_NAME(scale_fields_table_g)(d.data(), my_fields->metal_density,
                   &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2],
                   &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], imetal, &factor);

  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */

