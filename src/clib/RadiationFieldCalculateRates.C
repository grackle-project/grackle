/***********************************************************************
/
/  SETS THE MULTI-SPECIES RATES BASED ON THE EXTERNAL RADIATION FIELD
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  Pascal Paschos, Robert Harkness
/  date:       1st July 2002
/  modified2:  Pascal Paschos
/  date:       August, 2002	
/  modified3:  Robert Harkness - Killed 32-bit IBM C++ bug
/  date:       15 September 2002
/  modified4:  Pascal Paschos & Robert Harkness - parameter controls
/  date:       March, 2006
/  modified5:  Robert Harkness - removed C++ I/O
/  date:       February 29th, 2008
/  modified6:  Ji-hoon Kim - changes to work with non-cosmological runs
/  date:       November, 2009
/
/  PURPOSE:    
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

/* function prototypes */

int RadiationFieldCalculateRates(chemistry_data &my_chemistry,
                                 code_units &my_units, float a_value)
{
  /* Return if there is no radiation (rates should be all zero). */

  if (my_chemistry.RadiationFieldType == 0)
    return SUCCESS;

  /* Set units. */

  float Redshift = 1.0 / (a_value * my_units.a_units) - 1;
  if (!my_units.comoving_coordinates) {
    my_chemistry.RadiationRedshiftOn = Redshift+0.2;
    my_chemistry.RadiationRedshiftOff = 0.0;
    my_chemistry.RadiationRedshiftFullOn = Redshift+0.1;
    my_chemistry.RadiationRedshiftDropOff = 0.0;
  }
    
  double tbase1 = my_units.time_units;
  double xbase1 = my_units.length_units/(a_value * my_units.a_units);
  double dbase1 = my_units.density_units*POW(a_value * my_units.a_units, 3);
  double mh     = 1.67e-24;
  double CoolingUnits = (POW(my_units.a_units, 5) * xbase1*xbase1 * mh*mh) /
                        (POW(tbase1, 3) * dbase1);

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. 
     (Only used if my_chemistry.RadiationFieldType = 1 to 4 or equal to 12). */

  float Ramp = 0;

  if (Redshift < my_chemistry.RadiationRedshiftOn && 
      Redshift >= my_chemistry.RadiationRedshiftOff) {

    if (Redshift > my_chemistry.RadiationRedshiftFullOn)
      Ramp = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
	    (my_chemistry.RadiationRedshiftOn+my_chemistry.RadiationRedshiftFullOn)));
    else if (Redshift < my_chemistry.RadiationRedshiftDropOff)
      Ramp = (Redshift - my_chemistry.RadiationRedshiftOff + my_chemistry.f0to3*
	                 (my_chemistry.RadiationRedshiftDropOff - Redshift)) /
             (my_chemistry.RadiationRedshiftDropOff - 
	      my_chemistry.RadiationRedshiftOff);
    else
      Ramp = 1.0;

  }

  float RampX = 0.0;  // this is unused below
  float XRadRedShiftOn = 7.0;
  float XRadRedShiftOff = 0.0 ;
  float XRadRedShiftFullOn = 6.0 ;
  float XRadRedShiftDropOff = 0.0 ;
  float XRadRedShiftf0to3 = 0.1 ;
  
  if (Redshift < XRadRedShiftOn &&
      Redshift >= XRadRedShiftOff) {
    
    if (Redshift > XRadRedShiftFullOn)
      RampX = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
				   (XRadRedShiftOn + XRadRedShiftFullOn)));
    else if (Redshift < XRadRedShiftDropOff)
      RampX = (Redshift - XRadRedShiftOff + XRadRedShiftf0to3*
	       (XRadRedShiftDropOff - Redshift)) /
	(XRadRedShiftDropOff -
	 XRadRedShiftOff);
    else
      RampX = 1.0;
    
  }

  float exp_arg = -1.0 * POW(Redshift-2.3, 2);

  /* Above z = 2.3, increase the width of the Gaussian by a factor of
     3 for HI and HeI. */

  float beta2 = (my_chemistry.AdjustUVBackgroundHighRedshift && Redshift > 2.3) ? 9.0 : 1.0;

  /* ------------------------------------------------------------------ */
  /* 1) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5) */

  if (my_chemistry.RadiationFieldType == 1) {

    my_chemistry.k24 = 6.7e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                     * my_units.time_units * Ramp;
    my_chemistry.k25 = 6.3e-15 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     * my_units.time_units * Ramp;
    my_chemistry.k26 = 3.2e-13 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00/beta2) 
                     * my_units.time_units * Ramp;
    my_chemistry.piHI   = 4.7e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2) 
                     / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 8.2e-24 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00/beta2) 
                     / CoolingUnits * Ramp;
    my_chemistry.piHeII = 1.6e-25 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     / CoolingUnits * Ramp;
  }   
    
  /* ------------------------------------------------------------------ */
  /* 2) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.8) */

  if (my_chemistry.RadiationFieldType == 2) {
    my_chemistry.k24 = 5.6e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 * my_units.time_units * Ramp;
    my_chemistry.k25 = 3.2e-15 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.60)
                 * my_units.time_units * Ramp;
    my_chemistry.k26 = 4.8e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 * my_units.time_units * Ramp;
    my_chemistry.piHI   = 3.9e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 6.4e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/2.10/beta2)
                 / CoolingUnits * Ramp;
    my_chemistry.piHeII = 8.7e-26 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.70)
                 / CoolingUnits * Ramp;
  }

  /* ------------------------------------------------------------------ */
  /* 3) This used to be a modified version of (1) but with the HeII heating rate
     multiplied by 1.8. Since the parameter is now controlled externally this option
     is not needed altogether. Instead, we are going to use it as the 
     HM-12 background reconstructed and modified to match observations 
     by Kirkman & Tytler (2005). */ 
    

  if (my_chemistry.RadiationFieldType == 3) {
    my_chemistry.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * my_units.time_units * Ramp;
    my_chemistry.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;

    float KP_mod = 1.0 ; 

    if (Redshift >= 3 ) {
      KP_mod = 1.3 ;
    }
    if (Redshift > 2 && Redshift < 3) {
      KP_mod = (1.0+(Redshift-2.0)*0.3) ;
    }

    my_chemistry.k24 *= KP_mod ;
    my_chemistry.k25 *= KP_mod ;
    my_chemistry.k26 *= KP_mod ;
    my_chemistry.piHI *= KP_mod ; 
    my_chemistry.piHeI *= KP_mod ;
    my_chemistry.piHeII *= KP_mod ;

                                                                            
}

  /* ------------------------------------------------------------------ */
  /* 4) For the Haardt and Madau (2001) quasar spectrum (alpha_q = 1.57) 
        with X-ray Compton heating from Madau & Efstathiou (9902080). */

  if (my_chemistry.RadiationFieldType == 4) {

    my_chemistry.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * my_units.time_units * Ramp;
    my_chemistry.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;

    /* This is sigma_thompson * c * (effective <h \mu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \mu> is the
       average photon energy in keV, corrected for relativistic effects. 
       Eq.(4) and Eq.(11) of Madau & Efstathiou (1999) */

    float RedshiftXrayCutoff = 5.0;
    
    my_chemistry.comp_xray = 6.65e-25 * 3.0e10 * 
                        (31.8*POW(1.0+Redshift, 0.3333)/511.0) * 
                        (6.3e-5 * 1.6e-12) * 
                        POW(1.0 + Redshift, 4) * 
                        exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits; 
    /*
    my_chemistry.comp_xray = 6.65e-25 * 3.0e10 * 
                        (1.0/511.0e3) * 
                        (4.0 * 1.38e-16/1.6e-12) *
                        (6.3e-5 * 1.6e-12) * 
                        POW(1.0 + Redshift, 4) * 
                        exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits; */

    /* The effective temperature (in K).  Eq.(10) of Madau & Efstathiou (1999) 
       with U_xray(z=0) = 6.3e-5 and U_cmb(z=0) = 0.256 eV/cm3  */

    my_chemistry.temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
                         (4.0*1.38e-16) *
                         6.3e-5 * POW(1.0 + Redshift, 4) * 
                         exp(-POW(Redshift/RedshiftXrayCutoff, 2)) /
                         (0.256 * (1+Redshift));  

    /*
    my_chemistry.temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
       (4.0*1.38e-16); */  // <-- this seems to be wrong, Ji-hoon Kim
  }


  /* ------------------------------------------------------------------ */
  /* 5) Featureless power-law spectrum (see calc_rates.src). */

  if (my_chemistry.RadiationFieldType == 5)
    ;

  /* ------------------------------------------------------------------ */
  /* 8) An absorbed (hard) quasar-like spectrum plus molecular H constant
     photo-dissociation  (the rates are calculated in calc_rates.src). */

  if (my_chemistry.RadiationFieldType == 8) {

    /* Insert redshift-dependent term here */

    /* molecular hydrogen constant photo-dissociation
       rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands)
       Note: this is hard-coded to F_LW = 1e-21; flux normalization controls
       the hard-radiation component (i.e. > 13.6 eV)*/

    my_chemistry.k31 = 1.13e8 * 1.0e-21 * my_units.time_units;
  }

  /* ------------------------------------------------------------------ */
  /* 9) molecular hydrogen constant photo-dissociation only! 
     rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) */

  if (my_chemistry.RadiationFieldType == 9) {

      my_chemistry.k31 = 1.13e8 * my_chemistry.f3 * my_units.time_units;

  }

  /* ------------------------------------------------------------------ */
  /* 10 & 11) - internally-computed radiation field.  Most of the rates
     are calculated in RadiationFieldUpdate, but the comp_xray rate is
     calculated here, from the X-ray energy density and temperature. */

  // if (my_chemistry.RadiationFieldType >= 10 && my_chemistry.RadiationFieldType <= 11) {
  //   my_chemistry.comp_xray = 8.0/3.0*0.665e6/(9.1e0*3.0e0) *
  //     RadiationData.ComptonXrayEnergyDensity*1.38e2*
  //     (double(1.0e-30)/CoolingUnits);
  //   my_chemistry.temp_xray = RadiationData.ComptonXrayTemperature;
  // }

  /* ------------------------------------------------------------------ */
  /* 12) For the Haardt and Madau (2001) QSO+GAL (alpha_q = 1.57)        */

  if (my_chemistry.RadiationFieldType == 12) {

   
    my_chemistry.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * my_units.time_units * Ramp;
    my_chemistry.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * my_units.time_units * Ramp;
    my_chemistry.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    my_chemistry.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;
    }

 /* 13) Fits to the (unpublished) Haardt & Madau (2005) model that
        ships with Cloudy v8.00 (Michael Kuhlen, 06/02/2010).

        Like my_chemistry.RadiationFieldType=4, this model includes both a quasar
        and a galaxy contribution, but no X-ray Compton heating. */

  if (my_chemistry.RadiationFieldType == 13) {

    my_chemistry.k24 = 4.68e-12 * POW(1.0+Redshift, -0.592) * 
      exp( -0.7156 * POW(Redshift-2.292, 2.0) / 
	   (1.0 + 0.2009 * POW(Redshift+0.548, 2.0)) ) *
      my_units.time_units * Ramp;

    my_chemistry.k25 = 3.67e-14 * POW(1.0+Redshift, -0.258) * 
      exp( -1.1378 * POW(Redshift-1.875, 2.0) / 
	   (1.0 + 0.1727 * POW(Redshift-0.607, 2.0)) ) *
      my_units.time_units * Ramp;

    my_chemistry.k26 = 1.11e-11 * POW(1.0+Redshift, -1.451) * 
      exp( -0.7552 * POW(Redshift-2.954, 2.0) / 
	   (1.0 + 0.2625 * POW(Redshift+0.918, 2.0)) ) *
      my_units.time_units * Ramp;

    if (my_chemistry.primordial_chemistry > 1) {

      my_chemistry.k27 = 9.37e-10 * POW(1.0+Redshift, 2.526) *
	exp( 12.7560 * POW(Redshift+2.707, 2.0) /
	     (1.0 + -0.0301 * POW(Redshift+111.200, 2.0)) ) *
	my_units.time_units * Ramp;
      
      my_chemistry.k28 = 9.84e-12 * POW(1.0+Redshift, 2.115) * 
	exp( -1.9137 * POW(Redshift-1.748, 2.0) / 
	     (1.0 + 0.3896 * POW(Redshift+3.679, 2.0)) ) * 
	my_units.time_units * Ramp;
      
      my_chemistry.k29 = 6.87e-12 * POW(1.0+Redshift, -0.456) * 
	exp( -0.7601 * POW(Redshift-2.234, 2.0) / 
	     (1.0 + 0.1955 * POW(Redshift+0.655, 2.0)) ) * 
	my_units.time_units * Ramp;
      
      my_chemistry.k30 = 2.61e-12 * POW(1.0+Redshift, -2.198) * 
	exp( -0.5824 * POW(Redshift-3.811, 2.0) / 
	     (1.0 + 0.3241 * POW(Redshift+0.838, 2.0)) ) * 
	my_units.time_units * Ramp;
      
      /* LymanSawtoothSuppressionFactor is supposed to account for the
	 suppression of LW flux due to Lyman-series absorption (giving
	 a sawtooth pattern), a la Haiman & Abel, & Rees (2000).  This
	 is really only applicable when there is plenty of neutral
	 hydrogen is around, so I'm scaling it with Ramp. */
      float LymanSawtoothSuppressionFactor = 0.1 + 0.9 * Ramp;

      my_chemistry.k31 = 2.36e-12 * POW(1.0+Redshift, 1.185) * 
	exp( -0.2173 * POW(Redshift-3.211, 2.0) / 
	     (1.0 + 1.1852 * POW(Redshift+0.156, 2.0)) ) * 
	my_units.time_units * Ramp * LymanSawtoothSuppressionFactor;

    }

    my_chemistry.piHI = 4.09e-23 * POW(1.0+Redshift, -0.746) * 
      exp( -0.7469 * POW(Redshift-2.419, 2.0) / 
	   (1.0 + 0.2120 * POW(Redshift+0.686, 2.0)) )
      / CoolingUnits * Ramp;

    my_chemistry.piHeII = 1.28e-24 * POW(1.0+Redshift, -0.504) *
      exp( -1.0742 * POW(Redshift-1.889, 2.0) /
	   (1.0 + 0.1919 * POW(Redshift-0.518, 2.0)) )
      / CoolingUnits * Ramp;

    my_chemistry.piHeI = 4.86e-22 * POW(1.0+Redshift, -2.302) * 
      exp( -0.5250 * POW(Redshift-3.900, 2.0) / 
	   (1.0 + 0.3452 * POW(Redshift+0.673, 2.0)) )
      / CoolingUnits * Ramp;

 }

  /* ------------------------------------------------------------------ */
  /* 14) molecular hydrogen photo-dissociation only with a fit from
     the Wise & Abel (2005) model, updated for WMAP7, reionization at
     z=6.8.  Only valid for 6<z<50.  At z>50, set to tiny.  At z<6,
     set it to my_chemistry.f3.

     rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) */

  float tiny_number = 1.0e-20;
  if (my_chemistry.RadiationFieldType == 14) {
    float logJ;
    if (Redshift > 50.0)
      my_chemistry.k31 = tiny_number;
    else if (Redshift > 6.0) {
      logJ = -23.56688 + 4.56213e-1 * (1.0+Redshift) -
	2.67982e-2 * POW(1.0+Redshift, 2.0) + 
	5.88234e-4 * POW(1.0+Redshift, 3.0) -
	5.05576e-6 * POW(1.0+Redshift, 4.0);
      my_chemistry.k31 = 1.13e8 * POW(10.0,logJ) * my_units.time_units;  //*4.0*M_PI
    } else
      my_chemistry.k31 = 1.13e8 * my_chemistry.f3 * my_units.time_units;
  }

/* ------------------------------------------------------------------ */
  if (my_chemistry.RadiationFieldType < 0 || my_chemistry.RadiationFieldType > 14) {
    ENZO_VFAIL("my_chemistry.RadiationFieldType %"ISYM" not recognized.\n", 
	    my_chemistry.RadiationFieldType)
   }

  if (my_chemistry.AdjustUVBackground < 0 || my_chemistry.AdjustUVBackground > 2 ) {
   ENZO_VFAIL("my_chemistry.AdjustUVBackground Type %"ISYM" not recognized.\n",
            my_chemistry.AdjustUVBackground)
  }

  if (my_chemistry.AdjustUVBackground == 0 ) {
/* All rates are computed at their theoretical values */
    return SUCCESS;
  }

#ifdef NO_HEII_UVB
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("WARNING: Not using a radiation background for HeIII\n");
  my_chemistry.k25 = 0.0;
  my_chemistry.piHeII = 0.0;
#endif /* NO_HEII_UVB */

  if (my_chemistry.RadiationFieldType == 1 || my_chemistry.RadiationFieldType == 2 || 
      my_chemistry.RadiationFieldType == 3 || my_chemistry.RadiationFieldType == 4 || 
      my_chemistry.RadiationFieldType == 12 || my_chemistry.RadiationFieldType == 13) {
  
    if (my_chemistry.AdjustUVBackground == 1 ) {
      /* Increases the HM-type HeII photoheating rates by the default value of 1.8 */
      my_chemistry.SetUVBAmplitude = 1.0 ;
      my_chemistry.SetHeIIHeatingScale = 1.8 ;
      my_chemistry.piHeII *= (my_chemistry.SetUVBAmplitude * my_chemistry.SetHeIIHeatingScale);
    }
    if (my_chemistry.AdjustUVBackground == 2) {

      /* Alters the HM-type rates by the user supplied values.
	 If some are not supplied then it uses the default 
	 values defined in the previous if statement and set
	 in the SetDefaultGlobalValues.C  */

      my_chemistry.k24    *= my_chemistry.SetUVBAmplitude ;
      my_chemistry.k25    *= my_chemistry.SetUVBAmplitude ;
      my_chemistry.k26    *= my_chemistry.SetUVBAmplitude ;
      my_chemistry.piHI   *= my_chemistry.SetUVBAmplitude ;
      my_chemistry.piHeI  *= my_chemistry.SetUVBAmplitude ;
      my_chemistry.piHeII *= (my_chemistry.SetUVBAmplitude * my_chemistry.SetHeIIHeatingScale);
    }

  }

  return SUCCESS;
}
