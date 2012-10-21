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
#include "global_data.h"
#include "chemistry_data.h"
#include "code_units.h"

/* function prototypes */

int update_UVbackground_rates(chemistry_data &my_chemistry,
			      code_units &my_units, float a_value)
{
  /* Return if there is no radiation (rates should be all zero). */

  if (my_chemistry.UVbackground_type == 0)
    return SUCCESS;


  /* Check that UVbackground_type is allowed. */

  if (my_chemistry.UVbackground_type < -2 || my_chemistry.UVbackground_type> 4) {
    ENZO_VFAIL("my_chemistry.UVbackground_type %"ISYM" not recognized.\n", 
	    my_chemistry.UVbackground_type)
   }


  /* Set units. */

  float Redshift = 1.0 / (a_value * my_units.a_units) - 1;
  if (!my_units.comoving_coordinates) {
    my_chemistry.UVbackground_redshift_on = Redshift+0.2;
    my_chemistry.UVbackground_redshift_off = 0.0;
    my_chemistry.UVbackground_redshift_fullon = Redshift+0.1;
    my_chemistry.UVbackground_redshift_drop = 0.0;
  }
    
  double tbase1 = my_units.time_units;
  double xbase1 = my_units.length_units/(a_value * my_units.a_units);
  double dbase1 = my_units.density_units*POW(a_value * my_units.a_units, 3);
  double mh     = 1.67e-24;
  double CoolingUnits = (POW(my_units.a_units, 5) * xbase1*xbase1 * mh*mh) /
                        (POW(tbase1, 3) * dbase1);

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. */

  float Ramp = 0;

  if (Redshift < my_chemistry.UVbackground_redshift_on && 
      Redshift >= my_chemistry.UVbackground_redshift_off) {

    if (Redshift > my_chemistry.UVbackground_redshift_fullon)
      Ramp = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
	    (my_chemistry.UVbackground_redshift_on+my_chemistry.UVbackground_redshift_fullon)));
    else if (Redshift < my_chemistry.UVbackground_redshift_drop)
      Ramp = (Redshift - my_chemistry.UVbackground_redshift_off + 
	      (my_chemistry.UVbackground_redshift_drop - Redshift)) /
             (my_chemistry.UVbackground_redshift_drop - 
	      my_chemistry.UVbackground_redshift_off);
    else
      Ramp = 1.0;

  }


  /* ------------------------------------------------------------------ */
  /* Haardt & Madau (1996) quasar spectrum (alpha_q = 1.5) */
  /* *** DEPRECATED *** */

  if (my_chemistry.UVbackground_type == -1) {
    float exp_arg = -1.0 * POW(Redshift-2.3, 2);

    my_chemistry.k24 = 6.7e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                     * my_units.time_units * Ramp;
    my_chemistry.k25 = 6.3e-15 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     * my_units.time_units * Ramp;
    my_chemistry.k26 = 3.2e-13 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00) 
                     * my_units.time_units * Ramp;
    my_chemistry.piHI   = 4.7e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95) 
                     / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 8.2e-24 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00) 
                     / CoolingUnits * Ramp;
    my_chemistry.piHeII = 1.6e-25 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     / CoolingUnits * Ramp;
  }   
    
  /* ------------------------------------------------------------------ */
  /* Haardt & Madau (1996) quasar spectrum (alpha_q = 1.8) */
  /* *** DEPRECATED *** */

  if (my_chemistry.UVbackground_type == -2) {
    float exp_arg = -1.0 * POW(Redshift-2.3, 2);

    my_chemistry.k24 = 5.6e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * my_units.time_units * Ramp;
    my_chemistry.k25 = 3.2e-15 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.60)
                 * my_units.time_units * Ramp;
    my_chemistry.k26 = 4.8e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 * my_units.time_units * Ramp;
    my_chemistry.piHI   = 3.9e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95)
                 / CoolingUnits * Ramp;
    my_chemistry.piHeI  = 6.4e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/2.10)
                 / CoolingUnits * Ramp;
    my_chemistry.piHeII = 8.7e-26 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.70)
                 / CoolingUnits * Ramp;
  }

  /* ------------------------------------------------------------------ */
  /* Haardt & Madau (2001) quasar + galaxy (alpha_q = 1.57)        */

  if (my_chemistry.UVbackground_type == 1) {
   
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

  /* ------------------------------------------------------------------ */
  /* Haardt & Madau model (ca. 2005) that ships with Cloudy v.8.00 */
  /* [quasars and galaxies] */

  if (my_chemistry.UVbackground_type == 2) {

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
  /* Haardt & Madau (2012) [quasars and galaxies] */

  if (my_chemistry.UVbackground_type == 3) {

    // NOT YET

  }

  /* ------------------------------------------------------------------ */
  /* Faucher-Giguere et al. (2011) [quasars and galaxies] */

  if (my_chemistry.UVbackground_type == 4) {

    // NOT YET

  }



  /* Molecular hydrogen constant photo-dissociation */

  /* Note that k31 can be set by the UVbackground, in which case it is
     overwritten here if (LWbackground_intensity > 0.0). */

  if (my_chemistry.LWbackground_intensity > 0.0) 
    my_chemistry.k31 = 1.13e8 * my_chemistry.LWbackground_intensity * my_units.time_units;
  
  /* LWbackground_sawtooth_suppression is supposed to account for the
     suppression of LW flux due to Lyman-series absorption (giving a
     sawtooth pattern), a la Haiman & Abel, & Rees (2000).  This is
     really only applicable when there is plenty of neutral hydrogen
     is around, so I'm scaling it with Ramp: when Ramp==1 (fully
     re-ionized universe) there's no suppression. */

  if (my_chemistry.LWbackground_sawtooth_suppression) {
    float LymanSawtoothSuppressionFactor = 0.1 + 0.9 * Ramp;
  
    my_chemistry.k31 *= LymanSawtoothSuppressionFactor;
  }



  /* Compton X-ray heating */

  if (my_chemistry.Compton_xray_heating) {

    float RedshiftXrayCutoff = 5.0;

    /* This is sigma_thompson * c * (effective <h \nu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \nu>
       is the average photon energy in keV, corrected for relativistic
       effects.  Eq.(4) and Eq.(11) of Madau & Efstathiou (1999) */
    
    my_chemistry.comp_xray = 6.65e-25 * 3.0e10 * 
      (31.8*POW(1.0+Redshift, 0.3333)/511.0) * 
      (6.3e-5 * 1.6e-12) * 
      POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
      CoolingUnits; 
    
  
    /* The effective temperature (in K).  Eq.(10) of Madau &
       Efstathiou (1999) with U_xray(z=0) = 6.3e-5 and U_cmb(z=0) =
       0.256 eV/cm3 */
    
    my_chemistry.temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
      (4.0*1.38e-16) *
      6.3e-5 * POW(1.0 + Redshift, 4) * 
      exp(-POW(Redshift/RedshiftXrayCutoff, 2)) /
      (0.256 * (1+Redshift));  

  }



  return SUCCESS;
}
