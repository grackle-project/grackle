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

int calc_rates_dust_loc (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_C13 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_C20 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_C25 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_C30 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_F13 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_F15 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_F50 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_F80 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_P170(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_P200(int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);
int calc_rates_dust_Y19 (int iSN, chemistry_data *my_chemistry, chemistry_data_storage *my_rates, double kunit, double coolunit);


int calc_rates_dust(chemistry_data *my_chemistry,
                    chemistry_data_storage *my_rates,
                    code_units *my_units)
{

//-------kdMgSiO3  :    Mg   +    SiO + 2H2O -> MgSiO3  + 2H2I
//-------kdAC      :    C                    -> AC      
                                 
//-------kdSiM     :    Si                   -> SiM     
//-------kdFeM     :    Fe                   -> FeM     
//-------kdMg2SiO4 :   2Mg   +    SiO + 3H2O -> Mg2SiO4 + 3H2I
//-------kdFe3O4   :   3Fe   +   4H2O        -> Fe3O4   + 4H2I
//-------kdSiO2D   :    SiO2                 -> SiO2D   
//-------kdMgO     :    Mg   +    H2O        -> MgO     +  H2I
//-------kdFeS     :    Fe   +    S          -> FeS     
//-------kdAl2O3   :   2Al   +   3H2O        -> Al2O3   + 3H2I
                     
//-------kdreforg  : 0.5CO   + 0.5CH2 + 1.2N -> reforg (C:H:O:N = 1:1:0.5:1.2)
//-------kdvolorg  :    CO   +   2H2I        -> volorg (CH3OH)
//-------kdH2Oice  :    H2O                  -> H2Oice


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

      my_rates->gr_N    = calloc(2, sizeof(int));
      my_rates->grain_N = calloc(2, sizeof(int));

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
      int NSN, NTd, Nfd, Nmom;
      double Td0, fd0;
      double dTd;
      int iSN, iTd, imom, itab;

      NSN = 12;
      my_rates->SN0_N = NSN;

      my_rates->SN0_XC  = malloc(NSN * sizeof(double));
      my_rates->SN0_XO  = malloc(NSN * sizeof(double));
      my_rates->SN0_XMg = malloc(NSN * sizeof(double));
      my_rates->SN0_XAl = malloc(NSN * sizeof(double));
      my_rates->SN0_XSi = malloc(NSN * sizeof(double));
      my_rates->SN0_XS  = malloc(NSN * sizeof(double));
      my_rates->SN0_XFe = malloc(NSN * sizeof(double));

      my_rates->SN0_fC  = malloc(NSN * sizeof(double));
      my_rates->SN0_fO  = malloc(NSN * sizeof(double));
      my_rates->SN0_fMg = malloc(NSN * sizeof(double));
      my_rates->SN0_fAl = malloc(NSN * sizeof(double));
      my_rates->SN0_fSi = malloc(NSN * sizeof(double));
      my_rates->SN0_fS  = malloc(NSN * sizeof(double));
      my_rates->SN0_fFe = malloc(NSN * sizeof(double));

      my_rates->SN0_fSiM      = malloc(NSN * sizeof(double));
      my_rates->SN0_fFeM      = malloc(NSN * sizeof(double));
      my_rates->SN0_fMg2SiO4  = malloc(NSN * sizeof(double));
      my_rates->SN0_fMgSiO3   = malloc(NSN * sizeof(double));
      my_rates->SN0_fFe3O4    = malloc(NSN * sizeof(double));
      my_rates->SN0_fAC       = malloc(NSN * sizeof(double));
      my_rates->SN0_fSiO2D    = malloc(NSN * sizeof(double));
      my_rates->SN0_fMgO      = malloc(NSN * sizeof(double));
      my_rates->SN0_fFeS      = malloc(NSN * sizeof(double));
      my_rates->SN0_fAl2O3    = malloc(NSN * sizeof(double));
      my_rates->SN0_freforg   = malloc(NSN * sizeof(double));
      my_rates->SN0_fvolorg   = malloc(NSN * sizeof(double));
      my_rates->SN0_fH2Oice   = malloc(NSN * sizeof(double));

      for(iSN = 0; iSN < NSN; iSN++) {
        my_rates->SN0_XC [iSN] = 0.0;
        my_rates->SN0_XO [iSN] = 0.0;
        my_rates->SN0_XMg[iSN] = 0.0;
        my_rates->SN0_XAl[iSN] = 0.0;
        my_rates->SN0_XSi[iSN] = 0.0;
        my_rates->SN0_XS [iSN] = 0.0;
        my_rates->SN0_XFe[iSN] = 0.0;

        my_rates->SN0_fC [iSN] = 0.0;
        my_rates->SN0_fO [iSN] = 0.0;
        my_rates->SN0_fMg[iSN] = 0.0;
        my_rates->SN0_fAl[iSN] = 0.0;
        my_rates->SN0_fSi[iSN] = 0.0;
        my_rates->SN0_fS [iSN] = 0.0;
        my_rates->SN0_fFe[iSN] = 0.0;

        my_rates->SN0_fSiM     [iSN] = 0.0;
        my_rates->SN0_fFeM     [iSN] = 0.0;
        my_rates->SN0_fMg2SiO4 [iSN] = 0.0;
        my_rates->SN0_fMgSiO3  [iSN] = 0.0;
        my_rates->SN0_fFe3O4   [iSN] = 0.0;
        my_rates->SN0_fAC      [iSN] = 0.0;
        my_rates->SN0_fSiO2D   [iSN] = 0.0;
        my_rates->SN0_fMgO     [iSN] = 0.0;
        my_rates->SN0_fFeS     [iSN] = 0.0;
        my_rates->SN0_fAl2O3   [iSN] = 0.0;
        my_rates->SN0_freforg  [iSN] = 0.0;
        my_rates->SN0_fvolorg  [iSN] = 0.0;
        my_rates->SN0_fH2Oice  [iSN] = 0.0;
      }

      my_rates->SN0_r0SiM      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0FeM      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Mg2SiO4  = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0MgSiO3   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Fe3O4    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0AC       = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0SiO2D    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0MgO      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0FeS      = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0Al2O3    = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0reforg   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0volorg   = malloc(NSN * 3 * sizeof(double));
      my_rates->SN0_r0H2Oice   = malloc(NSN * 3 * sizeof(double));

      itab = 0;
      for(iSN = 0; iSN < NSN; iSN++) {
        for(imom = 0; imom < 3; imom++) {
          my_rates->SN0_r0SiM     [itab] = 0.0;
          my_rates->SN0_r0FeM     [itab] = 0.0;
          my_rates->SN0_r0Mg2SiO4 [itab] = 0.0;
          my_rates->SN0_r0MgSiO3  [itab] = 0.0;
          my_rates->SN0_r0Fe3O4   [itab] = 0.0;
          my_rates->SN0_r0AC      [itab] = 0.0;
          my_rates->SN0_r0SiO2D   [itab] = 0.0;
          my_rates->SN0_r0MgO     [itab] = 0.0;
          my_rates->SN0_r0FeS     [itab] = 0.0;
          my_rates->SN0_r0Al2O3   [itab] = 0.0;
          my_rates->SN0_r0reforg  [itab] = 0.0;
          my_rates->SN0_r0volorg  [itab] = 0.0;
          my_rates->SN0_r0H2Oice  [itab] = 0.0;
          itab++;
        }
      }

      NTd =            35;
      Td0 =     0.0000000;
      dTd =     0.1000000;
     Nmom =             4;

      my_rates->gr_Td = malloc(NTd * Nmom * sizeof(double));
  
      my_rates->gr_Size = NTd * Nmom;
      my_rates->gr_N[0] = Nmom;
      my_rates->gr_N[1] = NTd;
      my_rates->gr_dT   = dTd;
      for(iTd = 0; iTd < NTd; iTd++)
        my_rates->gr_Td[iTd] = Td0 + (double)iTd * dTd;

      my_rates->SN0_kpSiM      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeM      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMg2SiO4  = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgSiO3   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFe3O4    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAC       = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpSiO2D    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpMgO      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpFeS      = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpAl2O3    = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpreforg   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpvolorg   = malloc(NSN * Nmom * NTd * sizeof(double));
      my_rates->SN0_kpH2Oice   = malloc(NSN * Nmom * NTd * sizeof(double));

      itab = 0;
      for(iSN = 0; iSN < NSN; iSN++) {
        for(imom = 0; imom < Nmom; imom++) {
          for(iTd = 0; iTd < NTd; iTd++) {
            my_rates->SN0_kpSiM     [itab] = 0.0;
            my_rates->SN0_kpFeM     [itab] = 0.0;
            my_rates->SN0_kpMg2SiO4 [itab] = 0.0;
            my_rates->SN0_kpMgSiO3  [itab] = 0.0;
            my_rates->SN0_kpFe3O4   [itab] = 0.0;
            my_rates->SN0_kpAC      [itab] = 0.0;
            my_rates->SN0_kpSiO2D   [itab] = 0.0;
            my_rates->SN0_kpMgO     [itab] = 0.0;
            my_rates->SN0_kpFeS     [itab] = 0.0;
            my_rates->SN0_kpAl2O3   [itab] = 0.0;
            my_rates->SN0_kpreforg  [itab] = 0.0;
            my_rates->SN0_kpvolorg  [itab] = 0.0;
            my_rates->SN0_kpH2Oice  [itab] = 0.0;
            itab++;
          }
        }
      }

      ifunc = calc_rates_dust_loc ( 0, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_C13 ( 1, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_C20 ( 2, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_C25 ( 3, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_C30 ( 4, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_F13 ( 5, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_F15 ( 6, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_F50 ( 7, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_F80 ( 8, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_P170( 9, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_P200(10, my_chemistry, my_rates, kunit, coolunit);
      ifunc = calc_rates_dust_Y19 (11, my_chemistry, my_rates, kunit, coolunit);

  return SUCCESS;
}
