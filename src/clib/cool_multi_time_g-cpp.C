// See LICENSE file for license and copyright information

/// @file cool_multi_time_g-cpp.C
/// @brief Declares signature of cool_multi_time_g

// This file was initially generated automatically during conversion of the
// cool_multi_time_g function from FORTRAN to C++

#include <cstdio>
#include <vector>

#include "grackle.h"
#include "fortran_func_decls.h"
#include "internal_types.hpp"
#include "utils-cpp.hpp"

#include "cool_multi_time_g-cpp.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void cool_multi_time_g(
  gr_float* cooltime_data_, int* imetal, double* utem, double* uxyz,
  double* urho, chemistry_data* my_chemistry, chemistry_data_storage* my_rates,
  code_units* my_units, grackle_field_data* my_fields,
  photo_rate_storage* my_uvb_rates
)
{
  grackle::impl::View<gr_float***> cooltime(cooltime_data_, my_fields->grid_dimension[0], my_fields->grid_dimension[1], my_fields->grid_dimension[2]);

  // grain temperature
  std::vector<double> tSiM(my_fields->grid_dimension[0]);
  std::vector<double> tFeM(my_fields->grid_dimension[0]);
  std::vector<double> tMg2SiO4(my_fields->grid_dimension[0]);
  std::vector<double> tMgSiO3(my_fields->grid_dimension[0]);
  std::vector<double> tFe3O4(my_fields->grid_dimension[0]);
  std::vector<double> tAC(my_fields->grid_dimension[0]);
  std::vector<double> tSiO2D(my_fields->grid_dimension[0]);
  std::vector<double> tMgO(my_fields->grid_dimension[0]);
  std::vector<double> tFeS(my_fields->grid_dimension[0]);
  std::vector<double> tAl2O3(my_fields->grid_dimension[0]);
  std::vector<double> treforg(my_fields->grid_dimension[0]);
  std::vector<double> tvolorg(my_fields->grid_dimension[0]);
  std::vector<double> tH2Oice(my_fields->grid_dimension[0]);

  // Locals


  // Slice locals
 
  std::vector<double> p2d(my_fields->grid_dimension[0]);
  std::vector<double> tgas(my_fields->grid_dimension[0]);
  std::vector<double> mmw(my_fields->grid_dimension[0]);
  std::vector<double> tdust(my_fields->grid_dimension[0]);
  std::vector<double> metallicity(my_fields->grid_dimension[0]);
  std::vector<double> dust2gas(my_fields->grid_dimension[0]);
  std::vector<double> rhoH(my_fields->grid_dimension[0]);
  std::vector<double> edot(my_fields->grid_dimension[0]);


  // Iteration mask for multi_cool

  std::vector<gr_mask_type> itmask(my_fields->grid_dimension[0]);
  std::vector<gr_mask_type> itmask_metal(my_fields->grid_dimension[0]);

  // \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
  // =======================================================================
  const int dk = my_fields->grid_end[2] - my_fields->grid_start[2] + 1;
  const int dj = my_fields->grid_end[1] - my_fields->grid_start[1] + 1;

  // Convert densities from comoving to 'proper'

  if (my_units->comoving_coordinates == 1)  {

    gr_float factor = (gr_float)(std::pow(my_units->a_value,(-3)) );

     FORTRAN_NAME(scale_fields_g)(
      &my_chemistry->primordial_chemistry, imetal, &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry,
      &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->multi_metals, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation, &factor,
      &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2], &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
      my_fields->density, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
      my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
      my_fields->metal_density, my_fields->dust_density,
      my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
      my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
      my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
      my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
      my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
      my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
      my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
      my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
      my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
      my_fields->local_ISM_metal_density, my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
      my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
      my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);
 
  }

  // Loop over slices (in the k-direction)

  // parallelize the k and j loops with OpenMP
  // flat j and k loops for better parallelism
  //_// PORT: #ifdef _OPENMP
  //_// PORT: !$omp parallel do schedule(runtime) private(
  //_// PORT: !$omp&   p2d,
  //_// PORT: !$omp&   tgas,
  //_// PORT: !$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
  //_// PORT: !$omp&   edot,
  //_// PORT: !$omp&   itmask )
  //_// PORT: #endif
  //_// TODO_USE: OMP_PRAGMA("omp parallel")
  {
    //_// TODO: move remaining relevant declarations here to replace OMP private

    // each OMP thread separately initializes/allocates variables defined in
    // the current scope

    grackle::impl::LogTLinInterpScratchBuf logTlininterp_buf =
      grackle::impl::new_LogTLinInterpScratchBuf(my_fields->grid_dimension[0]);

    grackle::impl::Cool1DMultiScratchBuf cool1dmulti_buf =
      grackle::impl::new_Cool1DMultiScratchBuf(my_fields->grid_dimension[0]);

    // Cooling/heating slice locals
    grackle::impl::CoolHeatScratchBuf coolingheating_buf =
      grackle::impl::new_CoolHeatScratchBuf(my_fields->grid_dimension[0]);

    //_// TODO_USE: OMP_PRAGMA("omp for")
    for (int t = 0; t<=(dk * dj - 1); t++) {
      int k = t/dj      + my_fields->grid_start[2]+1;
      int j = grackle::impl::mod(t,dj) + my_fields->grid_start[1]+1;

      for (int i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        itmask[i-1] = MASK_TRUE;
      }

      // Compute the cooling rate
      int dummy_iter_arg=1;
      double comp1, comp2;

       FORTRAN_NAME(cool1d_multi_g)(
                   my_fields->density, my_fields->internal_energy, my_fields->x_velocity, my_fields->y_velocity, my_fields->z_velocity, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
                   &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2], &my_chemistry->NumberOfTemperatureBins,
                   &my_units->comoving_coordinates, &my_chemistry->primordial_chemistry, imetal, &my_chemistry->metal_cooling,
                   &my_chemistry->h2_on_dust, &my_chemistry->dust_chemistry, &my_chemistry->use_dust_density_field, &my_chemistry->dust_recombination_cooling,
                   &my_fields->grid_rank, &my_fields->grid_start[0], &my_fields->grid_end[0], &j, &k, &my_chemistry->ih2co, &my_chemistry->ipiht,
                   &dummy_iter_arg, &my_chemistry->photoelectric_heating,
                   &my_units->a_value, &my_chemistry->TemperatureStart, &my_chemistry->TemperatureEnd, &my_chemistry->SolarMetalFractionByMass, &my_chemistry->local_dust_to_gas_ratio,
                   utem, uxyz, &my_units->a_units, urho, &my_units->time_units,
                   &my_chemistry->Gamma, &my_chemistry->HydrogenFractionByMass,
                   my_rates->ceHI, my_rates->ceHeI, my_rates->ceHeII, my_rates->ciHI, my_rates->ciHeI,
                   my_rates->ciHeIS, my_rates->ciHeII, my_rates->reHII, my_rates->reHeII1,
                   my_rates->reHeII2, my_rates->reHeIII, my_rates->brem, &my_rates->comp, &my_rates->gammah,
                   &my_chemistry->interstellar_radiation_field, my_rates->regr, &my_rates->gamma_isrf, &my_uvb_rates->comp_xray, &my_uvb_rates->temp_xray,
                   &my_uvb_rates->piHI, &my_uvb_rates->piHeI, &my_uvb_rates->piHeII, &comp1, &comp2,
                   my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->metal_density, my_fields->dust_density,
                   my_rates->hyd01k, my_rates->h2k01, my_rates->vibh, my_rates->roth, my_rates->rotl,
                   coolingheating_buf.hyd01k, coolingheating_buf.h2k01, coolingheating_buf.vibh, coolingheating_buf.roth, coolingheating_buf.rotl,
                   my_rates->GP99LowDensityLimit, my_rates->GP99HighDensityLimit, coolingheating_buf.gpldl, coolingheating_buf.gphdl,
                   my_rates->HDlte, my_rates->HDlow, coolingheating_buf.hdlte, coolingheating_buf.hdlow,
                   my_rates->GAHI, my_rates->GAH2, my_rates->GAHe, my_rates->GAHp, my_rates->GAel,
                   my_rates->H2LTE, my_rates->gas_grain,
                   coolingheating_buf.ceHI, coolingheating_buf.ceHeI, coolingheating_buf.ceHeII, coolingheating_buf.ciHI, coolingheating_buf.ciHeI, coolingheating_buf.ciHeIS, coolingheating_buf.ciHeII,
                   coolingheating_buf.reHII, coolingheating_buf.reHeII1, coolingheating_buf.reHeII2, coolingheating_buf.reHeIII, coolingheating_buf.brem,
                   logTlininterp_buf.indixe, logTlininterp_buf.t1, logTlininterp_buf.t2, logTlininterp_buf.logtem, logTlininterp_buf.tdef, edot.data(),
                   tgas.data(), cool1dmulti_buf.tgasold, mmw.data(), p2d.data(), tdust.data(), metallicity.data(),
                   dust2gas.data(), rhoH.data(), cool1dmulti_buf.mynh, cool1dmulti_buf.myde,
                   cool1dmulti_buf.gammaha_eff, cool1dmulti_buf.gasgr_tdust, cool1dmulti_buf.regr,
                   &my_chemistry->self_shielding_method, &my_uvb_rates->crsHI, &my_uvb_rates->crsHeI, &my_uvb_rates->crsHeII,
                   &my_uvb_rates->k24, &my_uvb_rates->k26, &my_chemistry->use_radiative_transfer, my_fields->RT_heating_rate,
                   &my_chemistry->h2_optical_depth_approximation, &my_chemistry->cie_cooling, &my_chemistry->h2_cooling_rate, &my_chemistry->hd_cooling_rate, my_rates->cieco, coolingheating_buf.cieco,
                   &my_chemistry->cmb_temperature_floor, &my_chemistry->UVbackground, &my_chemistry->cloudy_electron_fraction_factor,
                   &my_rates->cloudy_primordial.grid_rank, my_rates->cloudy_primordial.grid_dimension,
                   my_rates->cloudy_primordial.grid_parameters[0], my_rates->cloudy_primordial.grid_parameters[1], my_rates->cloudy_primordial.grid_parameters[2], my_rates->cloudy_primordial.grid_parameters[3], my_rates->cloudy_primordial.grid_parameters[4],
                   &my_rates->cloudy_primordial.data_size, my_rates->cloudy_primordial.cooling_data, my_rates->cloudy_primordial.heating_data, my_rates->cloudy_primordial.mmw_data,
                   &my_rates->cloudy_metal.grid_rank, my_rates->cloudy_metal.grid_dimension,
                   my_rates->cloudy_metal.grid_parameters[0], my_rates->cloudy_metal.grid_parameters[1], my_rates->cloudy_metal.grid_parameters[2], my_rates->cloudy_metal.grid_parameters[3], my_rates->cloudy_metal.grid_parameters[4],
                   &my_rates->cloudy_metal.data_size, my_rates->cloudy_metal.cooling_data, my_rates->cloudy_metal.heating_data, &my_rates->cloudy_data_new,
                   &my_chemistry->use_volumetric_heating_rate, &my_chemistry->use_specific_heating_rate, my_fields->volumetric_heating_rate, my_fields->specific_heating_rate,
                   &my_chemistry->use_temperature_floor, &my_chemistry->temperature_floor_scalar, my_fields->temperature_floor,
                   &my_chemistry->use_isrf_field, my_fields->isrf_habing, itmask.data(), itmask_metal.data(),
                  &my_chemistry->metal_chemistry, &my_chemistry->grain_growth, &my_chemistry->use_primordial_continuum_opacity, &my_chemistry->tabulated_cooling_minimum_temperature,
                  my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
                  my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
                  my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
                  my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
                  my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
                  my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
                  my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
                  my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
                  my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
                  my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
                  my_rates->cieY06,
                  my_rates->LH2.props.dimension, &my_rates->LH2.props.data_size,
                  my_rates->LH2.props.parameters[0], my_rates->LH2.props.parameters[1], my_rates->LH2.props.parameters[2],
                  &my_rates->LH2.props.parameter_spacing[0], &my_rates->LH2.props.parameter_spacing[1], &my_rates->LH2.props.parameter_spacing[2], my_rates->LH2.data,
                  my_rates->LHD.props.dimension, &my_rates->LHD.props.data_size,
                  my_rates->LHD.props.parameters[0], my_rates->LHD.props.parameters[1], my_rates->LHD.props.parameters[2],
                  &my_rates->LHD.props.parameter_spacing[0], &my_rates->LHD.props.parameter_spacing[1], &my_rates->LHD.props.parameter_spacing[2], my_rates->LHD.data,
                  my_rates->LCI.props.dimension, &my_rates->LCI.props.data_size,
                  my_rates->LCI.props.parameters[0], my_rates->LCI.props.parameters[1], my_rates->LCI.props.parameters[2],
                  &my_rates->LCI.props.parameter_spacing[0], &my_rates->LCI.props.parameter_spacing[1], &my_rates->LCI.props.parameter_spacing[2], my_rates->LCI.data,
                  my_rates->LCII.props.dimension, &my_rates->LCII.props.data_size,
                  my_rates->LCII.props.parameters[0], my_rates->LCII.props.parameters[1], my_rates->LCII.props.parameters[2],
                  &my_rates->LCII.props.parameter_spacing[0], &my_rates->LCII.props.parameter_spacing[1], &my_rates->LCII.props.parameter_spacing[2], my_rates->LCII.data,
                  my_rates->LOI.props.dimension, &my_rates->LOI.props.data_size,
                  my_rates->LOI.props.parameters[0], my_rates->LOI.props.parameters[1], my_rates->LOI.props.parameters[2],
                  &my_rates->LOI.props.parameter_spacing[0], &my_rates->LOI.props.parameter_spacing[1], &my_rates->LOI.props.parameter_spacing[2], my_rates->LOI.data,
                  my_rates->LCO.props.dimension, &my_rates->LCO.props.data_size,
                  my_rates->LCO.props.parameters[0], my_rates->LCO.props.parameters[1], my_rates->LCO.props.parameters[2],
                  &my_rates->LCO.props.parameter_spacing[0], &my_rates->LCO.props.parameter_spacing[1], &my_rates->LCO.props.parameter_spacing[2], my_rates->LCO.data,
                  my_rates->LOH.props.dimension, &my_rates->LOH.props.data_size,
                  my_rates->LOH.props.parameters[0], my_rates->LOH.props.parameters[1], my_rates->LOH.props.parameters[2],
                  &my_rates->LOH.props.parameter_spacing[0], &my_rates->LOH.props.parameter_spacing[1], &my_rates->LOH.props.parameter_spacing[2], my_rates->LOH.data,
                  my_rates->LH2O.props.dimension, &my_rates->LH2O.props.data_size,
                  my_rates->LH2O.props.parameters[0], my_rates->LH2O.props.parameters[1], my_rates->LH2O.props.parameters[2],
                  &my_rates->LH2O.props.parameter_spacing[0], &my_rates->LH2O.props.parameter_spacing[1], &my_rates->LH2O.props.parameter_spacing[2], my_rates->LH2O.data,
                  my_rates->alphap.props.dimension, &my_rates->alphap.props.data_size,
                  my_rates->alphap.props.parameters[0], my_rates->alphap.props.parameters[1], &my_rates->alphap.props.parameter_spacing[0], &my_rates->alphap.props.parameter_spacing[1],
                  my_rates->alphap.data,
                  &my_chemistry->multi_metals, &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->use_multiple_dust_temperatures, &my_chemistry->dust_sublimation,
                  my_fields->local_ISM_metal_density,
                  my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
                  my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
                  my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density,
                  &my_rates->SN0_N,
                  my_rates->SN0_fSiM, my_rates->SN0_fFeM, my_rates->SN0_fMg2SiO4, my_rates->SN0_fMgSiO3,
                  my_rates->SN0_fFe3O4, my_rates->SN0_fAC, my_rates->SN0_fSiO2D, my_rates->SN0_fMgO,
                  my_rates->SN0_fFeS, my_rates->SN0_fAl2O3,
                  my_rates->SN0_freforg, my_rates->SN0_fvolorg, my_rates->SN0_fH2Oice,
                  my_rates->SN0_r0SiM, my_rates->SN0_r0FeM, my_rates->SN0_r0Mg2SiO4, my_rates->SN0_r0MgSiO3,
                  my_rates->SN0_r0Fe3O4, my_rates->SN0_r0AC, my_rates->SN0_r0SiO2D, my_rates->SN0_r0MgO,
                  my_rates->SN0_r0FeS, my_rates->SN0_r0Al2O3,
                  my_rates->SN0_r0reforg, my_rates->SN0_r0volorg, my_rates->SN0_r0H2Oice,
                  my_rates->gr_N, &my_rates->gr_Size, &my_rates->gr_dT, my_rates->gr_Td,
                  my_rates->SN0_kpSiM, my_rates->SN0_kpFeM, my_rates->SN0_kpMg2SiO4, my_rates->SN0_kpMgSiO3,
                  my_rates->SN0_kpFe3O4, my_rates->SN0_kpAC, my_rates->SN0_kpSiO2D, my_rates->SN0_kpMgO,
                  my_rates->SN0_kpFeS, my_rates->SN0_kpAl2O3,
                  my_rates->SN0_kpreforg, my_rates->SN0_kpvolorg, my_rates->SN0_kpH2Oice,
                  tSiM.data(), tFeM.data(), tMg2SiO4.data(), tMgSiO3.data(), tFe3O4.data(),
                  tAC.data(), tSiO2D.data(), tMgO.data(), tFeS.data(), tAl2O3.data(),
                  treforg.data(), tvolorg.data(), tH2Oice.data(),
                  my_rates->gas_grain2, &my_rates->gamma_isrf2
               );

      // Compute the cooling time on the slice
      //   (the gamma used here is the same as used to calculate the pressure
      //    in cool1d_multi_g)

      for (int i = my_fields->grid_start[0] + 1; i<=(my_fields->grid_end[0] + 1); i++) {
        double energy = std::fmax(p2d[i-1]/(my_chemistry->Gamma-1.), tiny_fortran_val);
        cooltime(i-1,j-1,k-1) = (gr_float)(energy/edot[i-1] );
      }

    }

    // cleanup temporaries
    grackle::impl::drop_LogTLinInterpScratchBuf(&logTlininterp_buf);
    grackle::impl::drop_Cool1DMultiScratchBuf(&cool1dmulti_buf);
    grackle::impl::drop_CoolHeatScratchBuf(&coolingheating_buf);

  }  // OMP_PRAGMA("omp parallel")

  // Convert densities back to comoving from 'proper'

  if (my_units->comoving_coordinates == 1)  {

    gr_float factor = (gr_float)(std::pow(my_units->a_value,3) );

     FORTRAN_NAME(scale_fields_g)(
      &my_chemistry->primordial_chemistry, imetal, &my_chemistry->use_dust_density_field, &my_chemistry->metal_chemistry,
      &my_chemistry->metal_abundances, &my_chemistry->dust_species, &my_chemistry->multi_metals, &my_chemistry->grain_growth, &my_chemistry->dust_sublimation, &factor,
      &my_fields->grid_start[0], &my_fields->grid_end[0], &my_fields->grid_start[1], &my_fields->grid_end[1], &my_fields->grid_start[2], &my_fields->grid_end[2], &my_fields->grid_dimension[0], &my_fields->grid_dimension[1], &my_fields->grid_dimension[2],
      my_fields->density, my_fields->e_density, my_fields->HI_density, my_fields->HII_density, my_fields->HeI_density, my_fields->HeII_density, my_fields->HeIII_density,
      my_fields->HM_density, my_fields->H2I_density, my_fields->H2II_density, my_fields->DI_density, my_fields->DII_density, my_fields->HDI_density, my_fields->DM_density, my_fields->HDII_density, my_fields->HeHII_density,
      my_fields->metal_density, my_fields->dust_density,
      my_fields->CI_density, my_fields->CII_density, my_fields->CO_density, my_fields->CO2_density,
      my_fields->OI_density, my_fields->OH_density, my_fields->H2O_density, my_fields->O2_density,
      my_fields->SiI_density, my_fields->SiOI_density, my_fields->SiO2I_density,
      my_fields->CH_density, my_fields->CH2_density, my_fields->COII_density, my_fields->OII_density,
      my_fields->OHII_density, my_fields->H2OII_density, my_fields->H3OII_density, my_fields->O2II_density,
      my_fields->Mg_density, my_fields->Al_density, my_fields->S_density, my_fields->Fe_density,
      my_fields->SiM_dust_density, my_fields->FeM_dust_density, my_fields->Mg2SiO4_dust_density, my_fields->MgSiO3_dust_density, my_fields->Fe3O4_dust_density,
      my_fields->AC_dust_density, my_fields->SiO2_dust_density, my_fields->MgO_dust_density, my_fields->FeS_dust_density, my_fields->Al2O3_dust_density,
      my_fields->ref_org_dust_density, my_fields->vol_org_dust_density, my_fields->H2O_ice_dust_density,
      my_fields->local_ISM_metal_density, my_fields->ccsn13_metal_density, my_fields->ccsn20_metal_density, my_fields->ccsn25_metal_density, my_fields->ccsn30_metal_density,
      my_fields->fsn13_metal_density, my_fields->fsn15_metal_density, my_fields->fsn50_metal_density, my_fields->fsn80_metal_density,
      my_fields->pisn170_metal_density, my_fields->pisn200_metal_density, my_fields->y19_metal_density);

  }

  return;
}

#ifdef __cplusplus
}  // extern "C"
#endif /* __cplusplus */
