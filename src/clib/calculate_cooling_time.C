/***********************************************************************
/
/  GRID CLASS (COMPUTE THE COOLING TIME FIELD)
/
/  written by: Greg Bryan
/  date:       April, 1995
/  modified1:  Elizabeth Harper-Clark, August 2009
/              added in CoolingModel parameter
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"
#include "phys_constants.h"
#include "fortran.def"
 
extern "C" void FORTRAN_NAME(cool_multi_time)(
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	gr_float *cooltime,
	gr_int *in, gr_int *jn, gr_int *kn, gr_int *nratec, gr_int *iexpand,
	hydro_method *imethod,
        gr_int *ispecies, gr_int *imetal, gr_int *imcool, gr_int *idust, gr_int *idim,
	gr_int *is, gr_int *js, gr_int *ks, gr_int *ie, gr_int *je, gr_int *ke, gr_int *ih2co,
	gr_int *ipiht, gr_int *igammah,
	gr_float *dt, gr_float *aye, gr_float *temstart, gr_float *temend,
	gr_float *utem, gr_float *uxyz, gr_float *uaye, gr_float *urho, gr_float *utim,
	gr_float *eta1, gr_float *eta2, gr_float *gamma, gr_float *z_solar,
	gr_float *ceHIa, gr_float *ceHeIa, gr_float *ceHeIIa, gr_float *ciHIa, gr_float *ciHeIa,
	gr_float *ciHeISa, gr_float *ciHeIIa, gr_float *reHIIa, gr_float *reHeII1a,
	gr_float *reHeII2a, gr_float *reHeIIIa, gr_float *brema, gr_float *compa, gr_float *gammaha,
	gr_float *comp_xraya, gr_float *comp_temp, gr_float *piHI, gr_float *piHeI, gr_float *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II, gr_float *DI, gr_float *DII, gr_float *HDI, gr_float *metal,
	gr_float *hyd01ka, gr_float *h2k01a, gr_float *vibha, gr_float *rotha, gr_float *rotla,
	gr_float *gpldl, gr_float *gphdl, gr_float *HDltea, gr_float *HDlowa,
	gr_float *gaHIa, gr_float *gaH2a, gr_float *gaHea, gr_float *gaHpa, gr_float *gaela,
	gr_float *gasgra, gr_float *metala, gr_int *n_xe, gr_float *xe_start, gr_float *xe_end,
	gr_float *inutot, gr_int *iradfield, gr_int *nfreq, gr_int *imetalregen,
	gr_int *iradshield, gr_float *avgsighp, gr_float *avgsighep, gr_float *avgsighe2p,
	gr_int *iradtrans, gr_float *photogamma,
	gr_int *ih2optical, gr_int *iciecool, gr_float *ciecoa,
 	gr_int *icmbTfloor, gr_int *iClHeat,
 	gr_float *clEleFra, gr_int *clGridRank, gr_int *clGridDim,
 	gr_float *clPar1, gr_float *clPar2, gr_float *clPar3, gr_float *clPar4, gr_float *clPar5,
 	gr_int *clDataSize, gr_float *clCooling, gr_float *clHeating);
 
int calculate_cooling_time(chemistry_data &my_chemistry,
			   code_units &my_units,
			   gr_float a_value, gr_float dt_value,
			   gr_int grid_rank, gr_int *grid_dimension,
			   gr_int *grid_start, gr_int *grid_end,
			   gr_float *density, gr_float *internal_energy,
			   gr_float *x_velocity, gr_float *y_velocity, gr_float  *z_velocity,
			   gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
			   gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
			   gr_float *H2I_density, gr_float *H2II_density,
			   gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
			   gr_float *e_density, gr_float *metal_density,
			   gr_float *cooling_time)
{
 
  /* Return if this doesn't concern us. */

  if (!my_chemistry.use_chemistry) return SUCCESS;

  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */

  // gr_float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection * 
  //                        double(LengthUnits) * CellWidth[0][0];
  // gr_float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection * 
  //                         double(LengthUnits) * CellWidth[0][0];
  // gr_float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection * 
  //                          double(LengthUnits) * CellWidth[0][0];
  gr_float HIShieldFactor, HeIShieldFactor, HeIIShieldFactor;
  HIShieldFactor = HeIShieldFactor = HeIIShieldFactor = 0.0;
  

  /* Call the fortran routine to solve cooling equations. */

  gr_int ierr = 0;
  gr_int i_method = 2;  // so total energy is internal energy
  gr_int MetalFieldPresent = 1;
  gr_float TemperatureUnits =  mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;
  gr_float DualEnergyFormalismEta1 = 0.0;
  gr_float DualEnergyFormalismEta2 = 0.0;
  gr_int RadiationFieldRecomputeMetalRates = 0;
  gr_int RadiativeTransfer = 0;
  gr_int RadiativeTransferCoupledRateSolver = 0;
  gr_int RTCoupledSolverIntermediateStep = 0;
  gr_int RadiativeTransferHydrogenOnly = 0;

  gr_float *kphHINum, *kphHeINum, *kphHeIINum, 
    *kdissH2INum, *gammaNum;

    FORTRAN_NAME(cool_multi_time)(
       density, internal_energy, x_velocity, y_velocity, z_velocity,
       e_density, HI_density, HII_density,
       HeI_density, HeII_density, HeIII_density,
       cooling_time,
       grid_dimension, grid_dimension+1, grid_dimension+2,
       &my_chemistry.NumberOfTemperatureBins, &my_units.comoving_coordinates,
       &i_method,
       &my_chemistry.primordial_chemistry, &MetalFieldPresent, &my_chemistry.metal_cooling, 
       &my_chemistry.h2_on_dust,
       &grid_rank, grid_start, grid_start+1, grid_start+2,
       grid_end, grid_end+1, grid_end+2,
       &my_chemistry.ih2co, &my_chemistry.ipiht, &my_chemistry.photoelectric_heating,
       &dt_value, &a_value, &my_chemistry.TemperatureStart,
       &my_chemistry.TemperatureEnd,
       &TemperatureUnits, &my_units.length_units, &my_units.a_units, &my_units.density_units, &my_units.time_units,
       &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &my_chemistry.Gamma,
       &my_chemistry.SolarMetalFractionByMass,
       my_chemistry.ceHI, my_chemistry.ceHeI, my_chemistry.ceHeII, my_chemistry.ciHI,
       my_chemistry.ciHeI,
       my_chemistry.ciHeIS, my_chemistry.ciHeII, my_chemistry.reHII,
       my_chemistry.reHeII1,
       my_chemistry.reHeII2, my_chemistry.reHeIII, my_chemistry.brem, &my_chemistry.comp, &my_chemistry.gammah,
       &my_chemistry.comp_xray, &my_chemistry.temp_xray,
       &my_chemistry.piHI, &my_chemistry.piHeI, &my_chemistry.piHeII,
       HM_density, H2I_density, H2II_density,
       DI_density, DII_density, HDI_density,
       metal_density,
       my_chemistry.hyd01k, my_chemistry.h2k01, my_chemistry.vibh,
       my_chemistry.roth, my_chemistry.rotl,
       my_chemistry.GP99LowDensityLimit, my_chemistry.GP99HighDensityLimit,
       my_chemistry.HDlte, my_chemistry.HDlow,
       my_chemistry.GAHI, my_chemistry.GAH2, my_chemistry.GAHe, my_chemistry.GAHp,
       my_chemistry.GAel, my_chemistry.gas_grain,
       my_chemistry.metals, &my_chemistry.NumberOfElectronFracBins, 
       &my_chemistry.ElectronFracStart, &my_chemistry.ElectronFracEnd,
       my_chemistry.Spectrum[0], &my_chemistry.RadiationFieldType,
       &my_chemistry.NumberOfFrequencyBins,
       &RadiationFieldRecomputeMetalRates,
       &my_chemistry.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
       &RadiativeTransfer, gammaNum, 
       &my_chemistry.h2_optical_depth_approximation, &my_chemistry.cie_cooling, my_chemistry.cieco,
       &my_chemistry.cmb_temperature_floor,
       &my_chemistry.include_metal_heating,
       &my_chemistry.CloudyElectronFractionFactor,
       &my_chemistry.CloudyCoolingGridRank,
       my_chemistry.CloudyCoolingGridDimension,
       my_chemistry.CloudyCoolingGridParameters[0],
       my_chemistry.CloudyCoolingGridParameters[1],
       my_chemistry.CloudyCoolingGridParameters[2],
       my_chemistry.CloudyCoolingGridParameters[3],
       my_chemistry.CloudyCoolingGridParameters[4],
       &my_chemistry.CloudyDataSize,
       my_chemistry.CloudyCooling, my_chemistry.CloudyHeating);
 
  return SUCCESS;
}
