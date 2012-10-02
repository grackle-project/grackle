/***********************************************************************
/
/  GRID CLASS (SOLVE THE COOLING/HEATING AND RATE EQUATIONS)
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  July, 2005 to solve cool and rate equations simultaneously
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"
#include "fortran.def"

/* function prototypes */

extern "C" void FORTRAN_NAME(solve_rate_cool)(
	float *d, float *e, float *ge, float *u, float *v, float *w, float *de,
	float *HI, float *HII, float *HeI, float *HeII, float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand, 
	hydro_method *imethod,
        int *idual, int *ispecies, int *imetal, int *imcool, int *idust, int *idim,
	int *is, int *js, int *ks, int *ie, int *je, int *ke, int *ih2co, 
	int *ipiht, int *igammah,
	float *dt, float *aye, float *temstart, float *temend,
	float *utem, float *uxyz, float *uaye, float *urho, float *utim,
	float *eta1, float *eta2, float *gamma, float *fh, float *dtoh,
	float *z_solar,
	float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, 
	float *k6a, float *k7a, float *k8a, float *k9a, float *k10a,
	float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a, 
	float *k15a,
        float *k16a, float *k17a, float *k18a, float *k19a, float *k22a,
	float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
	float *k30, float *k31,
	float *k50a, float *k51a, float *k52a, float *k53a, float *k54a,
	float *k55a, float *k56a,
	int *ndratec, float *dtemstart, float *dtemend, float *h2dusta, 
	float *ncrna, float *ncrd1a, float *ncrd2a,
	float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, 
	float *ciHeIa, 
	float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a, 
	float *reHeII2a, float *reHeIIIa, float *brema, float *compa, float *gammaha,
	float *comp_xraya, float *comp_temp, 
	float *piHI, float *piHeI, float *piHeII,
	float *HM, float *H2I, float *H2II, float *DI, float *DII, float *HDI,
	float *metal,
	float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
	float *gpldl, float *gphdl, float *HDltea, float *HDlowa,
	float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela,
	float *gasgra, float *metala, int *n_xe, float *xe_start, float *xe_end,
	float *inutot, int *iradtype, int *nfreq, int *imetalregen,
	int *iradshield, float *avgsighp, float *avgsighep, float *avgsighe2p,
	int *iradtrans, int *iradcoupled, int *iradstep, int *ierr,
	int *irt_honly,
	float *kphHI, float *kphHeI, float *kphHeII, 
	float *kdissH2I, float *photogamma,
	int *ih2optical, int *iciecool, int *ithreebody, float *ciecoa,
 	int *icmbTfloor, int *iClHeat,
 	float *clEleFra, int *clGridRank, int *clGridDim,
 	float *clPar1, float *clPar2, float *clPar3, float *clPar4, float *clPar5,
 	int *clDataSize, float *clCooling, float *clHeating);


int solve_chemistry(chemistry_data &my_chemistry,
                    code_units &my_units,
                    float a_value, float dt_value,
                    int grid_dimension, float *grid_size,
                    float *grid_start, float *grid_end,
                    float *density, float *internal_energy,
                    float *x_velocity, float *y_velocity, float  *z_velocity,
                    float *HI_density, float *HII_density, float *HM_density,
                    float *HeI_density, float *HeII_density, float *HeIII_density,
                    float *H2I_density, float *H2II_density,
                    float *DI_density, float *DII_density,
                    float *e_density, float *metal_density)
                      
{
  /* Return if this doesn't concern us. */
  if (!(MultiSpecies && RadiativeCooling)) return SUCCESS;

  /* Set up information for rates which depend on the radiation field. 
     Precompute factors for self shielding (this is the cross section * dx). */

  float HIShieldFactor = RadiationData.HIAveragePhotoHeatingCrossSection * 
                         double(LengthUnits) * CellWidth[0][0];
  float HeIShieldFactor = RadiationData.HeIAveragePhotoHeatingCrossSection * 
                          double(LengthUnits) * CellWidth[0][0];
  float HeIIShieldFactor = RadiationData.HeIIAveragePhotoHeatingCrossSection * 
                           double(LengthUnits) * CellWidth[0][0];

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  FORTRAN_NAME(solve_rate_cool)(
    density, internal_energy, internal_energy, x_velocity, y_velocity, z_velocity,
    e_density, HI_density, HII_density, 
    HeI_density, HeII_density, HeIII_density, 
    grid_dimension, grid_dimension+1, grid_dimension+2, 
    &my_chemistry.NumberOfTemperatureBins, &ComovingCoordinates, &HydroMethod, 
    &DualEnergyFormalism, &MultiSpecies, &MetalFieldPresent, &MetalCooling, 
    &H2FormationOnDust, 
    &GridRank, GridStartIndex, GridStartIndex+1, GridStartIndex+2, 
    GridEndIndex, GridEndIndex+1, GridEndIndex+2,
    &my_chemistry.ih2co, &my_chemistry.ipiht, &PhotoelectricHeating,
    &dtCool, &afloat, &my_chemistry.TemperatureStart, &my_chemistry.TemperatureEnd,
    &TemperatureUnits, &LengthUnits, &aUnits, &DensityUnits, &TimeUnits,
    &DualEnergyFormalismEta1, &DualEnergyFormalismEta2, &Gamma,
    &my_chemistry.HydrogenFractionByMass, &my_chemistry.DeuteriumToHydrogenRatio,
    &my_chemistry.SolarMetalFractionByMass,
    my_chemistry.k1, my_chemistry.k2, my_chemistry.k3, my_chemistry.k4, my_chemistry.k5, 
    my_chemistry.k6, my_chemistry.k7, my_chemistry.k8, my_chemistry.k9, my_chemistry.k10,
    my_chemistry.k11, my_chemistry.k12, my_chemistry.k13, my_chemistry.k13dd, my_chemistry.k14, 
    my_chemistry.k15, my_chemistry.k16,
    my_chemistry.k17, my_chemistry.k18, my_chemistry.k19, my_chemistry.k22,
    &my_chemistry.k24, &my_chemistry.k25, &my_chemistry.k26, &my_chemistry.k27,
    &my_chemistry.k28, &my_chemistry.k29, &my_chemistry.k30, &my_chemistry.k31,
    my_chemistry.k50, my_chemistry.k51, my_chemistry.k52, my_chemistry.k53,
    my_chemistry.k54, my_chemistry.k55, my_chemistry.k56,
    &my_chemistry.NumberOfDustTemperatureBins, &my_chemistry.DustTemperatureStart, 
    &my_chemistry.DustTemperatureEnd, my_chemistry.h2dust, 
    my_chemistry.n_cr_n, my_chemistry.n_cr_d1, my_chemistry.n_cr_d2,
    my_chemistry.ceHI, my_chemistry.ceHeI, my_chemistry.ceHeII, my_chemistry.ciHI,
    my_chemistry.ciHeI, 
    my_chemistry.ciHeIS, my_chemistry.ciHeII, my_chemistry.reHII, my_chemistry.reHeII1, 
    my_chemistry.reHeII2, my_chemistry.reHeIII, my_chemistry.brem, &my_chemistry.comp, &my_chemistry.gammah,
    &my_chemistry.comp_xray, &my_chemistry.temp_xray,
    &my_chemistry.piHI, &my_chemistry.piHeI, &my_chemistry.piHeII,
    HM_density, H2I_density, H2II_density,
    DI_density, DII_density, HDI_density,
    MetalPointer,
    my_chemistry.hyd01k, my_chemistry.h2k01, my_chemistry.vibh, my_chemistry.roth,my_chemistry.rotl,
    my_chemistry.GP99LowDensityLimit, my_chemistry.GP99HighDensityLimit, 
    my_chemistry.HDlte, my_chemistry.HDlow,
    my_chemistry.GAHI, my_chemistry.GAH2, my_chemistry.GAHe, my_chemistry.GAHp,
    my_chemistry.GAel, my_chemistry.gas_grain, 
    my_chemistry.metals, &my_chemistry.NumberOfElectronFracBins, 
    &my_chemistry.ElectronFracStart, &my_chemistry.ElectronFracEnd,
    RadiationData.Spectrum[0], &RadiationFieldType, 
    &RadiationData.NumberOfFrequencyBins, 
    &RadiationFieldRecomputeMetalRates,
    &RadiationData.RadiationShield, &HIShieldFactor, &HeIShieldFactor, &HeIIShieldFactor,
    &RadiativeTransfer, &RadiativeTransferCoupledRateSolver,
    &RTCoupledSolverIntermediateStep, &ierr,
    &RadiativeTransferHydrogenOnly,
    BaryonField[kphHINum], BaryonField[kphHeINum], BaryonField[kphHeIINum], 
    BaryonField[kdissH2INum], BaryonField[gammaNum],
    &H2OpticalDepthApproximation, &CIECooling, &ThreeBodyRate, my_chemistry.cieco,
    &CloudyCoolingData.CMBTemperatureFloor,
    &CloudyCoolingData.IncludeCloudyHeating,
    &CloudyCoolingData.CloudyElectronFractionFactor,
    &CloudyCoolingData.CloudyCoolingGridRank,
    CloudyCoolingData.CloudyCoolingGridDimension,
    CloudyCoolingData.CloudyCoolingGridParameters[0],
    CloudyCoolingData.CloudyCoolingGridParameters[1],
    CloudyCoolingData.CloudyCoolingGridParameters[2],
    CloudyCoolingData.CloudyCoolingGridParameters[3],
    CloudyCoolingData.CloudyCoolingGridParameters[4],
    &CloudyCoolingData.CloudyDataSize,
    CloudyCoolingData.CloudyCooling, CloudyCoolingData.CloudyHeating);

  return SUCCESS;

}
