#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define tiny 1.0e-20
#define huge 1.0e+20
#define tevk 1.1605e+4

extern int grackle_verbose;

int calc_rates_dust_loc(chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);


int calc_rates_dust(chemistry_data *my_chemistry,
                    chemistry_data_storage *my_rates,
                    code_units *my_units)
{

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

      int  ispecies    = my_chemistry->primordial_chemistry;
      int  igammah     = my_chemistry->photoelectric_heating;
      int  idust       = my_chemistry->h2_on_dust;
      int  idustall    = my_chemistry->dust_chemistry;
      int  nratec      = my_chemistry->NumberOfTemperatureBins;
      double  aye      = my_units->a_value;
      double  temstart = my_chemistry->TemperatureStart;
      double  temend   = my_chemistry->TemperatureEnd;
      int  casebrates  = my_chemistry->CaseBRecombination;
      int  threebody   = my_chemistry->three_body_rate;
      double  uxyz     = co_length_units;
      double  uaye     = my_units->a_units;
      double  urho     = co_density_units;
      double  utim     = my_units->time_units;

      int i,j,idt;
      double logttt, ttt, tev, logtev,
           xx, dum, tbase1, xbase1, kunit, coolunit,
           dbase1, dlogtem, kunit_3bdy, cierate,
           grain_coef, fgr, d_ttt, d_logttt, d_dlogtem,
           ttt2, ttt300, d_ttt2, tk9, lambdaHI, lambdaHeII,
           lambdaHeIII, grbeta;
      double tm, HDLR, HDLV, lt, t3, lt3;
      int anydust;
//
//    Set flag for dust-related options
//
      anydust = (idust > 0) || (idustall > 0);
//
//
// Get conversion units
//
//    t/x/dbase1 is the number (z dependant) that converts from the
//      dimensionless code units to physical units.  Also, in the
//      code aye = 1 at z=zinit, so to convert the usual a (=1 at z=0)
//      to a~ (written in the code as aye), we use a = a~*[a] 
//
      tbase1 = utim;
      xbase1 = uxyz/(aye*uaye);      // uxyz is [x]*a     
      dbase1 = urho*pow(aye*uaye, 3);// urho is [dens]/a^3
//
// 1) Set the dimensions of the (non-radiative) rate coefficients.  
//   Note that we have included the units that convert density to 
//   number density, so the rate equations should look like 
//   (in dimensionless units, hence the primes):
//
//      d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
//
//   where k~ is the dimenionless rate coefficients and d0-2~ are three
//    dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
//    dimensionless expansion coefficient (see above).
//
//   rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
//   rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
//   rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
//   so: k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] )  (~)
//   reminder: the number densities here are normalized with [dens] which
//             is not a constant (it has a factor a^3), so the number
//             densities must be converted from comoving to proper.
//
      kunit   = (pow(uaye, 3) * mh) / (dbase1 * tbase1);
      kunit_3bdy  = kunit * (pow(uaye, 3) * mh) / dbase1;
//
// 2) Set the dimension of the cooling coefficients (including constants)
//    (this equation has a rho because e is the specifi//energy, not
//     energy/unit volume).
//      delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
//      [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
//      delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
//    so L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
//      but [e] = ([a]*[x])**2 / [time]**2 and ([a] = 1 / (1 + zri) )
//     [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
//
      coolunit = (pow(uaye, 5) * pow(xbase1, 2) * pow(mh, 2)) / (pow(tbase1, 3) * dbase1);
//
//   Note: some of the coffiecients have only one power of n.  These
//         do not have the /a^3 factor, also they have units
//         [L1] = ([a]**2 * [x]**2 * mh) / [time]**3
//              = [L] * [dens] * [a]**3 / mh
//         This is done through the dom variable in cool.src
//        (some have three powers of n and they are different by the
//         reciprocal of the above factor multiplying [L]).
//
// 3) the units for the radiative rate coefficients is just 1/[time]
//
// 4) Energy transfer from gas to dust grains, following equation 2.15
//    of Hollenbach & McKee (1989).
//    Normalize to the HM89 dust to gas ratio.
      fgr = 0.009387;
      grain_coef = 1.2e-31 * pow(1.0e3, -0.5) / fgr;
//
// Compute log spacing in temperature
//
      ttt    = temstart;
      logttt = log(ttt);
      dlogtem= (log(temend) - log(temstart))/(double)(nratec-1);
//
// Compute log spacing in dust temperature
//
//    d_ttt     = dtemstart;
//    d_logttt  = log(d_ttt);
//    d_dlogtem = (log(dtemend) - log(dtemstart))/(double)(ndratec-1);

      if (ispecies == 0)
        return SUCCESS;

//    printf("%13.5e %13.5e %13.5e\n"
//     , tbase1
//     , xbase1
//     , dbase1
//    );
//
// Initialize constants to tiny
//
      int ifunc;
      int NTd, Nfd, Nmom;
      double Td0, fd0;
      double dTd;
      int iTd, ifd, itab;

      NTd =            35;
      Td0 =     0.0000000;
      dTd =     0.1000000;
     Nmom =             4;

      my_rates->gr_N  = malloc(2 * sizeof(int));
      my_rates->gr_Td = malloc(NTd * Nmom * sizeof(double));
  
      my_rates->gr_Size = NTd * Nmom;
      my_rates->gr_N[0] = Nmom;
      my_rates->gr_N[1] = NTd;
      my_rates->gr_dT   = dTd;
      for(iTd = 0; iTd < NTd; iTd++)
        my_rates->gr_Td[iTd] = Td0 + (double)iTd * dTd;

      ifunc = calc_rates_dust_loc(my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_C30(my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_F13(my_chemistry, my_rates, kunit, coolunit);

  return SUCCESS;
}
