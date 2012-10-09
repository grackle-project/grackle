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
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	float *cooltime,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
	hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co,
	int *ipiht, int *igammah,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma, float *z_solar,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
	float *comp_xraya, float *comp_temp, float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI, float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *gasgra, float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradfield, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, float *photogamma,
	int *ih2optical, int *iciecool, float *ciecoa,
 	int *icmbTfloor, int *iClHeat,
 	float *clEleFra, int *clGridRank, int *clGridDim,
 	float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
 	int *clDataSize, float *clCooling, float *clHeating);
 
int calculate_cooling_time(chemistry_data &my_chemistry,
			   code_units &my_units,
			   float a_value, float dt_value,
			   int grid_rank, int *grid_dimension,
			   int *grid_start, int *grid_end,
			   float *density, float *internal_energy,
			   float *x_velocity, float *y_velocity, float  *z_velocity,
			   float *HI_density, float *HII_density, float *HM_density,
			   float *HeI_density, float *HeII_density, float *HeIII_density,
			   float *H2I_density, float *H2II_density,
			   float *DI_density, float *DII_density, float *HDI_density,
			   float *e_density, float *metal_density,
			   float *cooling_time)
{
 
  /* Return if this doesn't concern us. */

  if (!my_chemistry.use_chemistry) return SUCCESS;

  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */

  // float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection * 
  //                        double(LengthUnits) * CellWidth[0][0];
  // float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection * 
  //                         double(LengthUnits) * CellWidth[0][0];
  // float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection * 
  //                          double(LengthUnits) * CellWidth[0][0];
  float HIShieldFactor, HeIShieldFactor, HeIIShieldFactor;
  HIShieldFactor = HeIShieldFactor = HeIIShieldFactor = 0.0;
  

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;
  int i_method = 2;  // so total energy is internal energy
  int i_dual = 1;    // doesn't matter
  int MetalFieldPresent = 1;
  float TemperatureUnits =  mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;
  float DualEnergyFormalismEta1 = 0.0;
  float DualEnergyFormalismEta2 = 0.0;
  int RadiationFieldRecomputeMetalRates = 0;
  int RadiativeTransfer = 0;
  int RadiativeTransferCoupledRateSolver = 0;
  int RTCoupledSolverIntermediateStep = 0;
  int RadiativeTransferHydrogenOnly = 0;

  float *kphHINum, *kphHeINum, *kphHeIINum, 
    *kdissH2INum, *gammaNum;

    FORTRAN_NAME(cool_multi_time)(
       density, internal_energy, internal_energy, x_velocity, y_velocity, z_velocity,
       e_density, HI_density, HII_density,
       HeI_density, HeII_density, HeIII_density,
       cooling_time,
       grid_dimension, grid_dimension+1, grid_dimension+2,
       &my_chemistry.NumberOfTemperatureBins, &my_units.comoving_coordinates,
       &i_method,
       &i_dual, &my_chemistry.primordial_chemistry, &MetalFieldPresent, &my_chemistry.metal_cooling, 
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
