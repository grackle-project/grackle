/***********************************************************************
/
/ Update UV background rates to current redshift
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              code_units *my_units)
{
  /* Return if there is no radiation (rates should be all zero). */

  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->primordial_chemistry == 0)
    return SUCCESS;

  /* Return if redshift is outside of table (rates should be all zero). */

  double Redshift = 1.0 / (my_units->a_value * my_units->a_units) - 1;
  if ( (Redshift < my_rates->UVbackground_table.zmin) ||
       (Redshift > my_rates->UVbackground_table.zmax) )
    return SUCCESS;

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. */

  double Ramp = 0;

  if (Redshift < my_chemistry->UVbackground_redshift_on && 
      Redshift >= my_chemistry->UVbackground_redshift_off) {

    if (Redshift > my_chemistry->UVbackground_redshift_fullon)
      Ramp = 0.5 - 0.5 * tanh(15.0 * (Redshift - 0.5 *
        (my_chemistry->UVbackground_redshift_on +
         my_chemistry->UVbackground_redshift_fullon)));
    else if (Redshift < my_chemistry->UVbackground_redshift_drop)
      Ramp = 0.5 - 0.5 * tanh(15.0 * (0.5 *
        (my_chemistry->UVbackground_redshift_drop +
         my_chemistry->UVbackground_redshift_off) - Redshift));
    else
      Ramp = 1.0;

  }


  /* Interpolate the UV background table. */
  
  // find interpolation index
  double *zvec = my_rates->UVbackground_table.z;
  double slope;
  int index=0;
  while (Redshift > zvec[index])
    index++;
  if(index == 0) index=1;
  if(index == my_rates->UVbackground_table.Nz) index--;

  // printf("index = %d, %.3f <= %.3f <= %.3f\n",index,zvec[index-1],Redshift,zvec[index]);

  // *** k24 ***
  slope = (my_rates->UVbackground_table.k24[index] -
           my_rates->UVbackground_table.k24[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->k24 = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.k24[index-1];

  // *** k25 ***
  slope = (my_rates->UVbackground_table.k25[index] -
           my_rates->UVbackground_table.k25[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->k25 = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.k25[index-1];

  // *** k26 ***
  slope = (my_rates->UVbackground_table.k26[index] -
           my_rates->UVbackground_table.k26[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->k26 = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.k26[index-1];

  if (my_chemistry->primordial_chemistry > 1) {

    // *** k27 ***
    slope = (my_rates->UVbackground_table.k27[index] -
             my_rates->UVbackground_table.k27[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->k27 = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.k27[index-1];

    // *** k28 ***
    slope = (my_rates->UVbackground_table.k28[index] -
             my_rates->UVbackground_table.k28[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->k28 = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.k28[index-1];

    // *** k29 ***
    slope = (my_rates->UVbackground_table.k29[index] -
             my_rates->UVbackground_table.k29[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->k29 = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.k29[index-1];

    // *** k30 ***
    slope = (my_rates->UVbackground_table.k30[index] -
             my_rates->UVbackground_table.k30[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->k30 = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.k30[index-1];

    // *** k31 ***
    slope = (my_rates->UVbackground_table.k31[index] -
             my_rates->UVbackground_table.k31[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->k31 = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.k31[index-1];

  }

  // *** piHI ***
  slope = (my_rates->UVbackground_table.piHI[index] -
           my_rates->UVbackground_table.piHI[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->piHI = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.piHI[index-1];

  // *** piHeII ***
  slope = (my_rates->UVbackground_table.piHeII[index] -
           my_rates->UVbackground_table.piHeII[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->piHeII = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.piHeII[index-1];

  // *** piHeI ***
  slope = (my_rates->UVbackground_table.piHeI[index] -
           my_rates->UVbackground_table.piHeI[index-1]) / (zvec[index] - zvec[index-1]);
  my_rates->piHeI = (Redshift - zvec[index-1]) * slope +
    my_rates->UVbackground_table.piHeI[index-1];

  //
  // Cross sections are read from table in cgs, and remain in cgs.
  // Analytic self-shielding method assumes cross sections in cgs units
  //

  // *** crsHI ***
  if (my_chemistry->self_shielding_method > 0){
    slope = (my_rates->UVbackground_table.crsHI[index] -
             my_rates->UVbackground_table.crsHI[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->hi_avg_crs = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.crsHI[index-1];

    // *** crsHeI ***
    slope = (my_rates->UVbackground_table.crsHeI[index] -
             my_rates->UVbackground_table.crsHeI[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->hei_avg_crs = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.crsHeI[index-1];

    // *** crsHeII ***
    slope = (my_rates->UVbackground_table.crsHeII[index] -
             my_rates->UVbackground_table.crsHeII[index-1]) / (zvec[index] - zvec[index-1]);
    my_rates->heii_avg_crs = (Redshift - zvec[index-1]) * slope +
      my_rates->UVbackground_table.crsHeII[index-1];
  }

  // Now convert the rates to code units.

  /* Get conversion units. */

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  double tbase1 = my_units->time_units;
  double xbase1 = co_length_units /
    (my_units->a_value * my_units->a_units);
  double dbase1 = co_density_units *
    POW(my_units->a_value * my_units->a_units, 3);
  double mh     = 1.67262171e-24;
  double ev2erg = 1.60217653e-12;
  /* compared to Enzo source, there's an additional factor of
     1/ev2erg here, because the heating rates are stored as eV/s. */
  double CoolingUnits = (POW(my_units->a_units, 5) * xbase1*xbase1 * mh*mh) /
    (POW(tbase1, 3) * dbase1) / ev2erg;


  my_rates->k24 *= my_units->time_units;
  my_rates->k25 *= my_units->time_units;
  my_rates->k26 *= my_units->time_units;
  
  if (my_chemistry->primordial_chemistry > 1) {
    my_rates->k27 *= my_units->time_units;
    my_rates->k28 *= my_units->time_units;
    my_rates->k29 *= my_units->time_units;
    my_rates->k30 *= my_units->time_units;
    my_rates->k31 *= my_units->time_units;
  }
    
  my_rates->piHI /= CoolingUnits;
  my_rates->piHeII /= CoolingUnits;
  my_rates->piHeI /= CoolingUnits;
  
  // Now apply the Ramp factor
  my_rates->k24 *= Ramp;
  my_rates->k25 *= Ramp;
  my_rates->k26 *= Ramp;
  if (my_chemistry->primordial_chemistry > 1) {
    my_rates->k27 *= Ramp;
    my_rates->k28 *= Ramp;
    my_rates->k29 *= Ramp;
    my_rates->k30 *= Ramp;
    my_rates->k31 *= Ramp;
  } 
  my_rates->piHI *= Ramp;
  my_rates->piHeII *= Ramp;
  my_rates->piHeI *= Ramp;

  /* Molecular hydrogen constant photo-dissociation */

  /* Note that k31 can be set by above by the UV background table, in
     which case it is overwritten here if (LWbackground_intensity >
     0.0). */

  if (my_chemistry->LWbackground_intensity > 0.0) 
    my_rates->k31 = 1.13e8 * my_chemistry->LWbackground_intensity * my_units->time_units;
  
  /* LWbackground_sawtooth_suppression is supposed to account for the
     suppression of LW flux due to Lyman-series absorption (giving a
     sawtooth pattern), a la Haiman & Abel, & Rees (2000).  This is
     really only applicable when there is plenty of neutral hydrogen
     is around, so I'm scaling it with Ramp: when Ramp==1 (fully
     re-ionized universe) there's no suppression. */

  if (my_chemistry->LWbackground_sawtooth_suppression) {
    double LymanSawtoothSuppressionFactor = 0.1 + 0.9 * Ramp;
  
    my_rates->k31 *= LymanSawtoothSuppressionFactor;
  }

  /* Compton X-ray heating */

  if (my_chemistry->Compton_xray_heating) {

    double RedshiftXrayCutoff = 5.0;

    /* This is sigma_thompson * c * (effective <h \nu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \nu>
       is the average photon energy in keV, corrected for relativistic
       effects.  Eq.(4) and Eq.(11) of Madau & Efstathiou (1999) */
    
    my_rates->comp_xray = 4.15e-13 * 3.0e10 *
      (31.8*POW(1.0+Redshift, 0.3333)/511.0) * 
      (6.3e-5 * 1.6e-12) * 
      POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
      CoolingUnits; 
    
  
    /* The effective temperature (in K).  Eq.(10) of Madau &
       Efstathiou (1999) with U_xray(z=0) = 6.3e-5 and U_cmb(z=0) =
       0.256 eV/cm3 */
    
    my_rates->temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
      (4.0*1.38e-16) *
      6.3e-5 * POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) /
      (0.256 * (1+Redshift));  

  }

  return SUCCESS;
}
