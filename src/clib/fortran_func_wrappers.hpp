// See LICENSE file for license and copyright information

/// @file fortran_func_wrappers.hpp
/// @brief Declares wrappers around routines that haven't been transcribed yet
///     from Fortran.
///
/// THIS FILE WILL BE DELETED ONCE WE COMPLETE TRANSCRIPTION
///
/// This file exists to aid the transcription process
/// - This file holds C++ functions that wrap untranscribed Fortran routines.
///   The idea is that these wrapper functions take reduced argument lists
///   (that may includes structs) before converting extended arg lists used
///   within the fortran subroutine.
/// - ideally the reduced argument list used by the wrapper function should
///   roughly approximate the final argument list that the routines will use
///   after transcription is complete.

#ifndef FORTRAN_FUNC_WRAPPERS_HPP
#define FORTRAN_FUNC_WRAPPERS_HPP

#ifndef __cplusplus
#error "This file must be read by a c++ compiler"
#endif


#include "dust/multi_grain_species/calc_grain_size_increment_1d.hpp"

#include "grackle.h"
#include "dust_props.hpp"
#include "fortran_func_decls.h"
#include "index_helper.h"
#include "inject_model/grain_metal_inject_pathways.hpp"
#include "internal_types.hpp"
#include "internal_units.h"
#include "LUT.hpp"
#include "opaque_storage.hpp"
#include "utils-cpp.hpp"

#include "step_rate_gauss_seidel.hpp"

// callers of these functions are generally expected to locally shorten the
// namespace name when they call these routines
namespace grackle::impl::fortran_wrapper {

inline void calc_all_tdust_gasgr_1d_g(
  double trad, double* tgas, double* tdust, double* metallicity,
  double* dust2gas, double* nh, double* gasgr_tdust, gr_mask_type* itmask_metal,
  double coolunit, double* gasgr, double* myisrf, double* kappa_tot,
  chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  grackle_field_data* my_fields, IndexRange idx_range,
  grackle::impl::GrainSpeciesCollection grain_temperatures,
  grackle::impl::GrainSpeciesCollection gas_grainsp_heatrate,
  grackle::impl::GrainSpeciesCollection grain_kappa,
  grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf,
  grackle::impl::InternalDustPropBuf internal_dust_prop_buf
) {

  // after transcription, we should obviously move this logic inside of
  // the transcribed function
  grackle::impl::GrainMetalInjectPathways* inject_pathway_props =
    my_rates->opaque_storage->inject_pathway_props;

  double dlog10Tdust = 0.0;
  double* log10Tdust_vals = nullptr;

  // NOTE: gr_N and gr_Size are historical names
  // -> they are pretty uninformative and should be changed!
  int gr_N[2] = {0, 0};
  int gr_Size = 0;
  if (inject_pathway_props != nullptr) {
    dlog10Tdust =
      inject_pathway_props->log10Tdust_interp_props.parameter_spacing[0];
    log10Tdust_vals =
      inject_pathway_props->log10Tdust_interp_props.parameters[0];

    gr_N[0] = inject_pathway_props->n_opac_poly_coef;
    gr_N[1] = static_cast<int>(
      inject_pathway_props->log10Tdust_interp_props.dimension[0]);
  };
  gr_Size = gr_N[0] * gr_N[1];


  FORTRAN_NAME(calc_all_tdust_gasgr_1d_g)(
    &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
    &my_chemistry->use_dust_density_field, &idx_range.i_start, &idx_range.i_end, &idx_range.jp1, &idx_range.kp1, &my_chemistry->local_dust_to_gas_ratio, &my_rates->gamma_isrf,
    &trad, my_rates->gas_grain, logTlininterp_buf.indixe, logTlininterp_buf.tdef, tgas, tdust,
    metallicity, dust2gas, nh, gasgr_tdust,
    itmask_metal,
    &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures,
    gr_N, &gr_Size, &dlog10Tdust, log10Tdust_vals,
    grain_temperatures.data[OnlyGrainSpLUT::SiM_dust], grain_temperatures.data[OnlyGrainSpLUT::FeM_dust], grain_temperatures.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_temperatures.data[OnlyGrainSpLUT::MgSiO3_dust], grain_temperatures.data[OnlyGrainSpLUT::Fe3O4_dust],
    grain_temperatures.data[OnlyGrainSpLUT::AC_dust], grain_temperatures.data[OnlyGrainSpLUT::SiO2_dust], grain_temperatures.data[OnlyGrainSpLUT::MgO_dust], grain_temperatures.data[OnlyGrainSpLUT::FeS_dust], grain_temperatures.data[OnlyGrainSpLUT::Al2O3_dust], grain_temperatures.data[OnlyGrainSpLUT::ref_org_dust],
    grain_temperatures.data[OnlyGrainSpLUT::vol_org_dust], grain_temperatures.data[OnlyGrainSpLUT::H2O_ice_dust], my_rates->gas_grain2, &my_rates->gamma_isrf2,
    &coolunit, gasgr, myisrf, internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Mg2SiO4_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::AC_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::FeS_dust],
    internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::Al2O3_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::ref_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_sigma_per_gas_mass.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.sigma_per_gas_mass_tot,
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeM_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::AC_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::SiO2_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::MgO_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::FeS_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::Al2O3_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::ref_org_dust],
    internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::vol_org_dust], internal_dust_prop_buf.grain_dyntab_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], internal_dust_prop_buf.dyntab_kappa_tot, grain_kappa.data[OnlyGrainSpLUT::SiM_dust], grain_kappa.data[OnlyGrainSpLUT::FeM_dust],
    grain_kappa.data[OnlyGrainSpLUT::Mg2SiO4_dust], grain_kappa.data[OnlyGrainSpLUT::MgSiO3_dust], grain_kappa.data[OnlyGrainSpLUT::Fe3O4_dust], grain_kappa.data[OnlyGrainSpLUT::AC_dust], grain_kappa.data[OnlyGrainSpLUT::SiO2_dust],
    grain_kappa.data[OnlyGrainSpLUT::MgO_dust], grain_kappa.data[OnlyGrainSpLUT::FeS_dust], grain_kappa.data[OnlyGrainSpLUT::Al2O3_dust], grain_kappa.data[OnlyGrainSpLUT::ref_org_dust], grain_kappa.data[OnlyGrainSpLUT::vol_org_dust],
    grain_kappa.data[OnlyGrainSpLUT::H2O_ice_dust], kappa_tot, gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiM_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeM_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Mg2SiO4_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgSiO3_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Fe3O4_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::AC_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::SiO2_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::MgO_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::FeS_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::Al2O3_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::ref_org_dust], gas_grainsp_heatrate.data[OnlyGrainSpLUT::vol_org_dust],
    gas_grainsp_heatrate.data[OnlyGrainSpLUT::H2O_ice_dust]
  );

}

inline void calc_temp1d_cloudy_g(
  double* rhoH, IndexRange idx_range, double* tgas, double* mmw, double dom,
  double zr, int imetal, cloudy_data cloudy_primordial, gr_mask_type* itmask,
  chemistry_data* my_chemistry, grackle_field_data* my_fields,
  InternalGrUnits internalu
) {
  FORTRAN_NAME(calc_temp1d_cloudy_g)(
    my_fields->density, my_fields->metal_density, my_fields->internal_energy, rhoH,
    &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_fields->grid_start[0], &my_fields->grid_end[0], &idx_range.jp1, &idx_range.kp1,
    tgas, mmw, &dom, &zr,
    &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd,
    &my_chemistry->Gamma, &internalu.utem, &imetal,
    &cloudy_primordial.grid_rank, cloudy_primordial.grid_dimension,
    cloudy_primordial.grid_parameters[0], cloudy_primordial.grid_parameters[1], cloudy_primordial.grid_parameters[2],
    &cloudy_primordial.data_size, cloudy_primordial.mmw_data,
    itmask
  );

}

/// Performs Gauss-Jordan elimination to solve the specified system of linear
/// equations and inverts the coefficient matrix.
///
/// In more detail, it solves the linear matrix `ax=b`, where `a` is the
/// square coefficient matrix, `b` is the right-hand side vector and `x`
/// is the solution vector
///
/// @param[in]     n The number of linear equations being solved
/// @param[in,out] coef_matrix An n by n column major array that initially
///    specifies the coefficient matrix. This function overwrites this matrix
///    with the inverted matrix.
/// @param[in,out] vec An n element array that initially specifies the
///    right-hand side vector. It's overwritten by the solution vector.
///
/// @retval 0 indicates success
/// @retval 1 indicates that the matrix is singular
///
/// > [!important]
/// > This appears to be taken directly from a routine provided in
/// > "Numerical Recipes" and slightly adapted. The original text allows
/// > `b` to hold multiple right-hand side vectors (and then we compute
/// > solutions for each right-hand side vector at the same time). In
/// > contrast, only allow 1 vector).
/// >
/// > The original Numerical Recipes snippet can be found here
/// > https://phys.uri.edu/nigh/NumRec/bookfpdf/f2-1.pdf
/// >
/// > Code from Numerical Recipes CANNOT be merged into Grackle (the
/// > licensing is fundamentally incompatible)
///
/// @todo
/// We need to replace this for reasons highlighted up above. When we do that,
/// we should account that the only place that calls this routine only needs
/// the solution to the system of equations (i.e. it does not need the
/// inverted matrix). For that reason, we may want to rename this to something
/// a little more generic. We should also consider information highlighted
/// within https://github.com/grackle-project/grackle/issues/255
inline int gaussj_g(int n, double* coef_matrix, double* vector) {
  int ierr;
  FORTRAN_NAME(gaussj_g)(&n, coef_matrix, vector, &ierr);
  return ierr;
}

/// wrapper for 1d interpolation
inline double interpolate_1d_g(
  double input1,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 1 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {
  double value;
  FORTRAN_NAME(interpolate_1d_g)(
    &input1, gridDim, gridPar1, &dgridPar1, &dataSize, dataField, &value
  );
  return value;
}

/// wrapper for 2d interpolation
inline double interpolate_2d_g(
  double input1, double input2,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 2 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  gr_i64 dataSize, const double* dataField
) {
  double value;
  FORTRAN_NAME(interpolate_2d_g)(
    &input1, &input2,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2,
    &dataSize, dataField, &value
  );
  return value;
}

/// Helper function used to implement interpolate_3dz_g
inline double interpolate_2df3d_g(
  double input1, double input3,
  const gr_i64* GRIMPL_RESTRICT gridDim, // 3 elements
  const double* GRIMPL_RESTRICT gridPar1, double dgridPar1,
  gr_i64 index2,
  const double* GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize,
  const double* GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_2df3d_g)(
    &input1, &input3,
    gridDim,
    gridPar1, &dgridPar1, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// Similar to interpolate_3d_g except index2 is calculated ahead of time
/// because it corresponds to redshift, which will not change for the entire
/// grid during the current calculation.
inline double interpolate_3dz_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, gr_i64 index2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField,
  gr_i64 end_int
) {

  double value;
  FORTRAN_NAME(interpolate_3dz_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &index2, gridPar3, &dgridPar3,
    &dataSize, dataField, &end_int, &value
  );
  return value;

}

/// wraps 3d interpolator functions
inline double interpolate_3d_g(
  double input1, double input2, double input3,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 3 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_3d_g)(
    &input1, &input2, &input3,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 4d interpolator functions
inline double interpolate_4d_g(
  double input1, double input2, double input3, double input4,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 4 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_4d_g)(
    &input1, &input2, &input3, &input4,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4,
    &dataSize, dataField, &value
  );
  return value;

}

/// wraps 5d interpolation routine
inline double interpolate_5d_g(
  double input1, double input2, double input3, double input4, double input5,
  const gr_i64 * GRIMPL_RESTRICT gridDim, // 5 elements
  const double * GRIMPL_RESTRICT gridPar1, double dgridPar1,
  const double * GRIMPL_RESTRICT gridPar2, double dgridPar2,
  const double * GRIMPL_RESTRICT gridPar3, double dgridPar3,
  const double * GRIMPL_RESTRICT gridPar4, double dgridPar4,
  const double * GRIMPL_RESTRICT gridPar5, double dgridPar5,
  gr_i64 dataSize, const double * GRIMPL_RESTRICT dataField
) {

  double value;
  FORTRAN_NAME(interpolate_5d_g)(
    &input1, &input2, &input3, &input4, &input5,
    gridDim,
    gridPar1, &dgridPar1, gridPar2, &dgridPar2, gridPar3, &dgridPar3,
    gridPar4, &dgridPar4, gridPar5, &dgridPar5,
    &dataSize, dataField, &value
  );
  return value;

}

} // namespace grackle::impl::fortran_wrapper

#endif /* FORTRAN_FUNC_WRAPPERS_HPP */

