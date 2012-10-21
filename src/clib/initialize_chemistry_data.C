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
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h" 
#include "phys_constants.h"

int initialize_cloudy_data(chemistry_data &my_chemistry,
                           code_units &my_units, gr_float a_value);

extern "C" void FORTRAN_NAME(calc_rates)(
     gr_int *nratec, gr_float *aye, gr_float *temstart, gr_float *temend, 
     gr_int *iradtype, gr_int *casebrates, gr_int *threebody,
     gr_float *utem, gr_float *uxyz, gr_float *uaye, gr_float *urho, gr_float *utim,
     gr_float *ceHIa, gr_float *ceHeIa, gr_float *ceHeIIa, gr_float *ciHIa, gr_float *ciHeIa,
     gr_float *ciHeISa, gr_float *ciHeIIa, gr_float *reHIIa, gr_float *reHeII1a,
     gr_float *reHeII2a, gr_float *reHeIIIa, gr_float *brema, gr_float *compa, 
     gr_float *gammahacgs, gr_float *gammaha,
     gr_float *hyd01ka, gr_float *h2k01a, gr_float *vibha, gr_float *rotha, gr_float *rotla,
     gr_float *gpldl, gr_float *gphdl, gr_float *hdlte, gr_float *hdlow, gr_float *hdcool, gr_float *cieco,
     gr_float *gaHIa, gr_float *gaH2a, gr_float *gaHea, gr_float *gaHpa, gr_float *gaela, gr_float *gasgr, 
     gr_float *k1a, gr_float *k2a, gr_float *k3a, gr_float *k4a, gr_float *k5a, gr_float *k6a,
        gr_float *k7a, gr_float *k8a, gr_float *k9a, gr_float *k10a,
     gr_float *k11a, gr_float *k12a, gr_float *k13a, gr_float *k13dda, gr_float *k14a,
        gr_float *k15a, gr_float *k16a, gr_float *k17a, gr_float *k18a,
     gr_float *k19a, gr_float *k20a, gr_float *k21a, gr_float *k22, gr_float *k23,
     gr_float *k50, gr_float *k51, gr_float *k52, gr_float *k53, gr_float *k54, gr_float *k55,
        gr_float *k56, gr_int *ndratec, gr_float *dtemstart, gr_float *dtemend, gr_float *h2dusta, 
     gr_float *ncrca, gr_float *ncrd1a, gr_float *ncrd2a, gr_int *ioutput);


// character strings
EXTERN char outfilename[];

 
int initialize_chemistry_data(chemistry_data &my_chemistry,
                              code_units &my_units, gr_float a_value)
{

  fprintf(stderr, "Initializing chemistry data.\n");

  /* Allocate CoolData space for rates. */
 
  my_chemistry.ceHI    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ceHeI   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ceHeII  = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHI    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeI   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeIS  = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.ciHeII  = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHII   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeII1 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeII2 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.reHeIII = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.brem    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.hyd01k  = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.h2k01   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.vibh    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.roth    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.rotl    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GP99LowDensityLimit  = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GP99HighDensityLimit = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDlte   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDlow   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.HDcool  = new gr_float[my_chemistry.NumberOfTemperatureBins*5];
  my_chemistry.cieco   = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHI    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAH2    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHe    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAHp    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.GAel    = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.gas_grain = new gr_float[my_chemistry.NumberOfTemperatureBins];

  /* Allocate space in my_chemistry for rates. */
 
  my_chemistry.k1 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k2 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k3 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k4 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k5 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k6 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k7 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k8 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k9 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k10 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k11 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k12 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k13 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k13dd = new gr_float[my_chemistry.NumberOfTemperatureBins*7];
  my_chemistry.k14 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k15 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k16 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k17 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k18 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k19 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k20 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k21 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k22 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k23 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k50 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k51 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k52 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k53 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k54 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k55 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.k56 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.h2dust = new gr_float[my_chemistry.NumberOfTemperatureBins * 
			      my_chemistry.NumberOfDustTemperatureBins];
  my_chemistry.n_cr_n = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.n_cr_d1 = new gr_float[my_chemistry.NumberOfTemperatureBins];
  my_chemistry.n_cr_d2 = new gr_float[my_chemistry.NumberOfTemperatureBins]; 

  gr_int ioutput = 1;
  gr_float temperature_units = mh*POW(my_units.length_units/
                                   my_units.time_units,2)/kboltz;

  /* Call FORTRAN routine to do the hard work. */
 
  FORTRAN_NAME(calc_rates)(
     &my_chemistry.NumberOfTemperatureBins, &a_value, &my_chemistry.TemperatureStart,
        &my_chemistry.TemperatureEnd,
        &my_chemistry.UVbackground_type, &my_chemistry.CaseBRecombination, &my_chemistry.three_body_rate,
     &temperature_units, &my_units.length_units, &my_units.a_units, 
     &my_units.density_units, &my_units.time_units,
     my_chemistry.ceHI, my_chemistry.ceHeI, my_chemistry.ceHeII, my_chemistry.ciHI,
        my_chemistry.ciHeI,
     my_chemistry.ciHeIS, my_chemistry.ciHeII, my_chemistry.reHII,
        my_chemistry.reHeII1,
     my_chemistry.reHeII2, my_chemistry.reHeIII, my_chemistry.brem, &my_chemistry.comp, 
     &my_chemistry.photoelectric_heating_rate, &my_chemistry.gammah,
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
