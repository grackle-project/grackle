/***********************************************************************
/
/  INITIALIZE THE MULTI-SPECIES RATES
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  Dan Reynolds, July 2010; added case-B recombination rates
/  modified2:  Britton Smith, October 2010; moved reading/writing of 
/              parameters to Read/WriteParameterFile.
/
/  PURPOSE:
/    For multi-species runs (with cooling), initialize both the
/      CoolData and RateData rate tables.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

int initialize_cloudy_data(chemistry_data &my_chemistry,
                           code_units &my_units, float a_value);

extern "C" void FORTRAN_NAME(calc_rates)(
     int *nratec, float *aye, float *temstart, float *temend, float *alpha0,
     float *f3, int *iradtype, int *casebrates, int *threebody,
     float *utem, float *uxyz, float *uaye, float *urho, float *utim,
     float *ceHIa, float *ceHeIa, float *ceHeIIa, float *ciHIa, float *ciHeIa,
     float *ciHeISa, float *ciHeIIa, float *reHIIa, float *reHeII1a,
     float *reHeII2a, float *reHeIIIa, float *brema, float *compa, 
     float *gammahacgs, float *gammaha,
     float *piHI, float *piHeI, float *piHeII,
     float *hyd01ka, float *h2k01a, float *vibha, float *rotha, float *rotla,
     float *gpldl, float *gphdl, float *hdlte, float *hdlow, float *hdcool, float *cieco,
     float *gaHIa, float *gaH2a, float *gaHea, float *gaHpa, float *gaela, float *gasgr, 
     float *k1a, float *k2a, float *k3a, float *k4a, float *k5a, float *k6a,
        float *k7a, float *k8a, float *k9a, float *k10a,
     float *k11a, float *k12a, float *k13a, float *k13dda, float *k14a,
        float *k15a, float *k16a, float *k17a, float *k18a,
     float *k19a, float *k20a, float *k21a, float *k22, float *k23,
     float *k24, float *k25, float *k26, float *k27, float *k28, float *k29,
        float *k30, float *k31,
     float *k50, float *k51, float *k52, float *k53, float *k54, float *k55,
        float *k56, int *ndratec, float *dtemstart, float *dtemend, float *h2dusta, 
     float *ncrca, float *ncrd1a, float *ncrd2a, int *ioutput);


// character strings
EXTERN char outfilename[];

 
int initialize_chemistry_data(chemistry_data &my_chemistry,
                              code_units &my_units, float a_value)
{

  fprintf(stderr, "Initializing chemistry data.\n");

  /* Allocate CoolData space for rates. */
 
  my_chemistry.ceHI    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ceHeI   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ceHeII  = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHI    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeI   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeIS  = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeII  = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHII   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeII1 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeII2 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeIII = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.brem    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.hyd01k  = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.h2k01   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.vibh    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.roth    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.rotl    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GP99LowDensityLimit  = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GP99HighDensityLimit = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDlte   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDlow   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDcool  = new float[my_chemistry.NumberOfTemperatureBins*5];
  my_chemistry.cieco   = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHI    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAH2    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHe    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHp    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAel    = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.gas_grain = new float[my_chemistry.NumberOfTemperatureBins];

  /* Allocate space in my_chemistry for rates. */
 
  my_chemistry.k1 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k2 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k3 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k4 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k5 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k6 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k7 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k8 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k9 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k10 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k11 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k12 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k13 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k13dd = new float[my_chemistry.NumberOfTemperatureBins*7];
  my_chemistry.k14 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k15 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k16 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k17 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k18 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k19 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k20 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k21 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k22 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k23 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k50 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k51 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k52 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k53 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k54 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k55 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k56 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.h2dust = new float[my_chemistry.NumberOfTemperatureBins * 
			      my_chemistry.NumberOfDustTemperatureBins];
  my_chemistry.n_cr_n = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.n_cr_d1 = new float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.n_cr_d2 = new float[my_chemistry.NumberOfTemperatureBins]; 

  int ioutput = 1;
  float temperature_units = mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;

  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates)(
     &my_chemistry.NumberOfTemperatureBins, &a_value, &my_chemistry.TemperatureStart,
        &my_chemistry.TemperatureEnd, &my_chemistry.alpha0, &my_chemistry.f3,
        &my_chemistry.RadiationFieldType, &my_chemistry.CaseBRecombination, &my_chemistry.three_body_rate,
     &temperature_units, &my_units.length_units, &my_units.a_units, 
     &my_units.density_units, &my_units.time_units,
     my_chemistry.ceHI, my_chemistry.ceHeI, my_chemistry.ceHeII, my_chemistry.ciHI,
        my_chemistry.ciHeI,
     my_chemistry.ciHeIS, my_chemistry.ciHeII, my_chemistry.reHII,
        my_chemistry.reHeII1,
     my_chemistry.reHeII2, my_chemistry.reHeIII, my_chemistry.brem, &my_chemistry.comp, 
     &my_chemistry.photoelectric_heating_rate, &my_chemistry.gammah,
     &my_chemistry.piHI, &my_chemistry.piHeI, &my_chemistry.piHeII,
     my_chemistry.hyd01k, my_chemistry.h2k01, my_chemistry.vibh, my_chemistry.roth,
        my_chemistry.rotl,
     my_chemistry.GP99LowDensityLimit, my_chemistry.GP99HighDensityLimit,
        my_chemistry.HDlte, my_chemistry.HDlow, my_chemistry.HDcool, my_chemistry.cieco,
     my_chemistry.GAHI, my_chemistry.GAH2, my_chemistry.GAHe, my_chemistry.GAHp,
        my_chemistry.GAel, my_chemistry.gas_grain, 
     my_chemistry.k1, my_chemistry.k2, my_chemistry.k3, my_chemistry.k4, my_chemistry.k5,
        my_chemistry.k6, my_chemistry.k7, my_chemistry.k8, my_chemistry.k9, my_chemistry.k10,
     my_chemistry.k11, my_chemistry.k12, my_chemistry.k13, my_chemistry.k13dd, my_chemistry.k14,
        my_chemistry.k15, my_chemistry.k16, my_chemistry.k17, my_chemistry.k18,
     my_chemistry.k19, my_chemistry.k20, my_chemistry.k21, my_chemistry.k22, my_chemistry.k23,
     &my_chemistry.k24, &my_chemistry.k25, &my_chemistry.k26, &my_chemistry.k27, &my_chemistry.k28,
        &my_chemistry.k29, &my_chemistry.k30, &my_chemistry.k31,
     my_chemistry.k50, my_chemistry.k51, my_chemistry.k52, my_chemistry.k53, my_chemistry.k54,
        my_chemistry.k55, my_chemistry.k56, 
     &my_chemistry.NumberOfDustTemperatureBins, &my_chemistry.DustTemperatureStart, 
     &my_chemistry.DustTemperatureEnd, my_chemistry.h2dust, 
     my_chemistry.n_cr_n, my_chemistry.n_cr_d1, my_chemistry.n_cr_d2, &ioutput);

  /* Initialize Cloudy cooling, even if not being used. */
  /* If not used, this will just initialize some data structues. */
  if (initialize_cloudy_data(my_chemistry, my_units, a_value) == FAIL) {
    fprintf(stderr, "Error in initialize_cloudy_data.");
    return FAIL;
  }

  return SUCCESS;
}
