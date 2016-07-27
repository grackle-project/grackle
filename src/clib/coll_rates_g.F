c=======================================================================
c/////////////////////  SUBROUTINE COLL_RATES_G  \\\\\\\\\\\\\\\\\\\\\\\
c
      SUBROUTINE coll_rates_g(T, k1, k2, k3, k4, k5, k6, k7, k8, k9,
     $     k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k23,
     $     kunit, casebrates)
c
c  COMPUTE MULTISPECIES COLLISIONAL RATES
c
c  written by: Tom Abel
c  date:       
c  modified1: Feb, 2000 by Greg Bryan; adapted to AMR
c  modified2: July, 2010 by Dan Reynolds; added case-B recombination rates
c
c  PURPOSE:
c    Computes various collisional rates (from Tom Abels web page)
c
c  UNITS:
c    cgs / kunit (kunit is a normalizing factor)
c
c  PARAMETERS:
c
c  INPUTS:
C     T is the gas temperature in Kelvin
c     kunit is a normalizing factor that (i.e. outputted rates are
c           cgs / kunit).
c
c  OUTPUTS:
c     k1-k19: rates as given below
c
C     the coefficient numbering is as in Abel etal 1997, NewA, 2.
C ---1:--       HI    + e   -> HII   + 2e
C ---2:--       HII   + e   -> H     + p
C ---3:--       HeI   + e   -> HeII  + 2e
C ---4:--       HeII  + e   -> HeI   + p
C ---5:--       HeII  + e   -> HeIII + 2e
C ---6:--       HeIII + e   -> HeII  + p
C ---7:--       HI    + e   -> HM    + p
C ---8:--       HM    + HI  -> H2I*  + e
C ---9:--       HI    + HII -> H2II  + p
C ---10--       H2II  + HI  -> H2I*  + HII
C ---11--       H2I   + HII -> H2II  + H
C ---12--       H2I   + e   -> 2HI   + e
C ---13--       H2I   + H   -> 3H
C ---14--       HM    + e   -> HI    + 2e
C ---15--       HM    + HI  -> 2H    + e
C ---16--       HM    + HII -> 2HI
C ---17--       HM    + HII -> H2II  + e
C ---18--       H2II  + e   -> 2HI
C ---19--       H2II  + HM  -> HI    + H2I
c
c-----------------------------------------------------------------------
c
#define USE_SAVIN2004
c
      implicit NONE
#include "grackle_fortran_types.def"
c
c  Arguments
c
      real*8 k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13,
     $       k14, k15, k16, k17, k18, k19, k23
      real*8 kunit, T
      integer casebrates
c
c  Parameters
c
c
c  Locals
c
      real*8 log_T, log_T_eV, T_ev, tk9
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
c=======================================================================
c
      
C     ------- Compute various values of T.
      log_T = LOG(T)
      T_eV = T/11605._DKIND
      log_T_eV = log(T_eV)
      
      k1 = exp(-32.71396786375_DKIND 
     &     + 13.53655609057_DKIND*log_T_eV
     &     - 5.739328757388_DKIND*log_T_eV**2 
     &     + 1.563154982022_DKIND*log_T_eV**3
     &     - 0.2877056004391_DKIND*log_T_eV**4
     &     + 0.03482559773736999_DKIND*log_T_eV**5
     &     - 0.00263197617559_DKIND*log_T_eV**6
     &     + 0.0001119543953861_DKIND*log_T_eV**7
     &     - 2.039149852002e-6_DKIND*log_T_eV**8) / kunit
      
      IF (T_eV .GT. 0.8) THEN
         k3 = exp(-44.09864886561001_DKIND
     &        + 23.91596563469_DKIND*log_T_eV
     &        - 10.75323019821_DKIND*log_T_eV**2
     &        + 3.058038757198_DKIND*log_T_eV**3
     &        - 0.5685118909884001_DKIND*log_T_eV**4
     &        + 0.06795391233790001_DKIND*log_T_eV**5
     &        - 0.005009056101857001_DKIND*log_T_eV**6
     &        + 0.0002067236157507_DKIND*log_T_eV**7
     &        - 3.649161410833e-6_DKIND*log_T_eV**8) / kunit
         
         k4 = (1.54e-9_DKIND*(1._DKIND+0.3_DKIND / 
     &        exp(8.099328789667_DKIND/T_eV))
     &        / (exp(40.49664394833662_DKIND/T_eV)*T_eV**1.5_DKIND)
     &        + 3.92e-13_DKIND/T_eV**0.6353_DKIND) / kunit
         
         k5 = exp(-68.71040990212001_DKIND
     &        + 43.93347632635_DKIND*log_T_eV
     &        - 18.48066993568_DKIND*log_T_eV**2
     &        + 4.701626486759002_DKIND*log_T_eV**3
     &        - 0.7692466334492_DKIND*log_T_eV**4
     &        + 0.08113042097303_DKIND*log_T_eV**5
     &        - 0.005324020628287001_DKIND*log_T_eV**6
     &        + 0.0001975705312221_DKIND*log_T_eV**7
     &        - 3.165581065665e-6_DKIND*log_T_eV**8) / kunit
      ELSE
         k1 = max(tiny, k1)
         k3 = tiny
         k4 = 3.92e-13_DKIND/T_eV**0.6353_DKIND / kunit
         k5 = tiny
      ENDIF

c     redefine k4 if case B recombination rates are requested
      if (casebrates.eq.1) then
         k4 = 1.26e-14_DKIND * (5.7067e5_DKIND/T)**(0.75_DKIND) / kunit
      endif

c     set HII recombination rate as either case A or case B
      IF (casebrates.eq.1) THEN
         IF (T .lt. 1.0d9) then
            k2 = 4.881357e-6_DKIND*T**(-1.5_DKIND) 
     &           * (1._DKIND+1.14813e2_DKIND
     &           * T**(-0.407_DKIND))**(-2.242_DKIND) / kunit
         ELSE
            k2 = tiny
         ENDIF
      ELSE
         IF ( T .GT. 5500._DKIND ) THEN
            k2 = exp(-28.61303380689232_DKIND
     &           - 0.7241125657826851_DKIND*log_T_eV
     &           - 0.02026044731984691_DKIND*log_T_eV**2
     &           - 0.002380861877349834_DKIND*log_T_eV**3
     &           - 0.0003212605213188796_DKIND*log_T_eV**4
     &           - 0.00001421502914054107_DKIND*log_T_eV**5
     &           + 4.989108920299513e-6_DKIND*log_T_eV**6
     &           + 5.755614137575758e-7_DKIND*log_T_eV**7
     &           - 1.856767039775261e-8_DKIND*log_T_eV**8
     &           - 3.071135243196595e-9_DKIND*log_T_eV**9) / kunit
         ELSE
            k2 = k4
         ENDIF
      ENDIF
      
c     set HeIII recombination rate as either case A or case B
      IF (casebrates.eq.1) THEN
         IF (T .lt. 1.0e9_DKIND) then
            k6 = 7.8155e-5_DKIND*T**(-1.5_DKIND) 
     &           * (1._DKIND+2.0189e2_DKIND
     &           * T**(-0.407_DKIND))**(-2.242_DKIND) / kunit
         ELSE
            k6 = tiny
         ENDIF
      ELSE
         k6 = 3.36e-10_DKIND/sqrt(T)/(T/1.e3_DKIND)**0.2_DKIND
     &        / (1._DKIND+(T/1.e6_DKIND)**0.7_DKIND) / kunit
      ENDIF
c
c     H- formation by radiative association: HI + e -> H2I + gamma.
c     Fit from Stancil, Lepp & Dalgarno (1998, ApJ, 509, 1), based
c     on photodetachment cross-section of Wishart (1979, MNRAS, 187, P59)       
c      
      k7 = 3.0e-16_DKIND * (T/3.e2_DKIND)**0.95_DKIND 
     &   * exp(-T/9.32e3_DKIND) / kunit
c
c     H- destruction by associative detachment with HI: HM + HI -> H2I + e
c     This fit is from Kreckel et al (2010, Science, 329, 69), based on their
c     experimental measurements of the process. 
c     
      k8 = 1.35e-9_DKIND*(T**9.8493e-2_DKIND + 3.2852e-1_DKIND
     &   * T**5.5610e-1_DKIND + 2.771e-7_DKIND * T**2.1826_DKIND)
     &   / (1._DKIND + 6.191e-3_DKIND * T**1.0461_DKIND
     &   + 8.9712e-11_DKIND * T**3.0424_DKIND
     &   + 3.2576e-14_DKIND * T**3.7741_DKIND) / kunit
      
c     H2+ formation by radiative association: HI + HII -> H2II + gamma
c     This fit is from Latif et al (2015, MNRAS, 446, 3163) and is valid
c     for temperatures in the range 1 < T < 32000 K. At higher temperatures,
c     we keep the value of the rate coefficient constant; note that its 
c     behaviour in this regime is unlikely to be important
c
      IF (T .GT. 3.2e4_DKIND) THEN
        tk9 = 3.2e4_DKIND
      ELSE
        tk9 = T
      ENDIF
      k9 = 1d1**(-18.20_DKIND  - 3.194_DKIND * log10(tk9)
     &          + 1.786_DKIND  * log10(tk9)**2
     &          - 0.2072_DKIND * log10(tk9)**3 ) / kunit
      
      k10 = 6.0e-10_DKIND / kunit
      
      IF (T_eV .GT. 0.3_DKIND) THEN
         k13 = 1.0670825e-10_DKIND*T_eV**2.012_DKIND/
     &        (exp(4.463_DKIND/T_eV)*(1._DKIND+0.2472_DKIND
     &        * T_eV)**3.512_DKIND) / kunit
         
#ifdef USE_SAVIN2004
         k11 = (exp(-21237.15_DKIND/T) * 
     $       (- 3.3232183e-07_DKIND
     $        + 3.3735382e-07_DKIND * log_T
     $        - 1.4491368e-07_DKIND * log_T**2
     $        + 3.4172805e-08_DKIND * log_T**3
     $        - 4.7813720e-09_DKIND * log_T**4
     $        + 3.9731542e-10_DKIND * log_T**5
     $        - 1.8171411e-11_DKIND * log_T**6
     $        + 3.5311932e-13_DKIND * log_T**7)) / kunit
#else /* Abel et al. (1997) */
         k11 = exp(-24.24914687731536_DKIND
     &        + 3.400824447095291_DKIND*log_T_eV
     &        - 3.898003964650152_DKIND*log_T_eV**2
     &        + 2.045587822403071_DKIND*log_T_eV**3
     &        - 0.5416182856220388_DKIND*log_T_eV**4
     &        + 0.0841077503763412_DKIND*log_T_eV**5
     &        - 0.007879026154483455_DKIND*log_T_eV**6
     &        + 0.0004138398421504563_DKIND*log_T_eV**7
     &        - 9.36345888928611e-6_DKIND*log_T_eV**8) / kunit
#endif /* USE_SAVIN2004 */
         

c Revised rate coefficient for collisional dissociation of H2 by electrons,
c from Trevisan & Tennyson (2002, Plasma Phys. Cont. Fus., 44, 1263). 
c Assumes all of the H2 molecules are in the v=0 vibrational ground state
      k12 = 4.4886e-9_DKIND*T**0.109127_DKIND*exp(-101858._DKIND/T) /
     &        kunit
      ELSE
         k13 = tiny 
         k11 = tiny
         k12 = tiny
      ENDIF
      
      IF (T_eV .GT. 0.04_DKIND) THEN
         k14 = exp(-18.01849334273_DKIND
     &        + 2.360852208681_DKIND*log_T_eV
     &        - 0.2827443061704_DKIND*log_T_eV**2
     &        + 0.01623316639567_DKIND*log_T_eV**3
     &        - 0.03365012031362999_DKIND*log_T_eV**4
     &        + 0.01178329782711_DKIND*log_T_eV**5
     &        - 0.001656194699504_DKIND*log_T_eV**6
     &        + 0.0001068275202678_DKIND*log_T_eV**7
     &        - 2.631285809207e-6_DKIND*log_T_eV**8) / kunit
      ELSE
         k14 =  tiny
      ENDIF
      
      IF (T_eV .GT. 0.1_DKIND) THEN
         k15 = exp(-20.37260896533324_DKIND
     &        + 1.139449335841631_DKIND*log_T_eV
     &        - 0.1421013521554148_DKIND*log_T_eV**2
     &        + 0.00846445538663_DKIND*log_T_eV**3
     &        - 0.0014327641212992_DKIND*log_T_eV**4
     &        + 0.0002012250284791_DKIND*log_T_eV**5
     &        + 0.0000866396324309_DKIND*log_T_eV**6
     &        - 0.00002585009680264_DKIND*log_T_eV**7
     &        + 2.4555011970392e-6_DKIND*log_T_eV**8
     &        - 8.06838246118e-8_DKIND*log_T_eV**9) / kunit
      ELSE
         k15 = 2.56e-9_DKIND*T_eV**1.78186_DKIND / kunit
      ENDIF
c
c     H- destruction by mutual neutralization with HII: HM + HII -> 2H
c     Fit from Croft et al (1999, MNRAS, 304, 327), based on the cross-section
c     of Fussen & Kubach (1986, J. Phys. B, 18 L31). This fit agrees well
c     with recent measurements of the rate by X. Urbain & collaborators, to
c     within experimental error
c      
      k16 = 2.4e-6_DKIND*(1._DKIND+T/2e4_DKIND)/sqrt(T) / kunit
      
      k17 = 1.e-8_DKIND*T**(-0.4_DKIND) / kunit
      IF (T .GT. 1.0e4_DKIND)
     &     k17=4.0e-4_DKIND*T**(-1.4_DKIND)*exp(-15100._DKIND/T) / kunit
      
      k18 = 1.e-8_DKIND / kunit
      If (T .GT. 617._DKIND)
     $     k18 = 1.32e-6_DKIND * T**(-0.76_DKIND) / kunit
      k19 = 5.e-7_DKIND*sqrt(100._DKIND/T) / kunit
      k23 = ((8.125e-8_DKIND / sqrt(T))
     $     * exp(-52000._DKIND/T)
     $     * (1._DKIND - exp(-6000._DKIND/T))) / kunit
      k23 = max(k23, tiny)
      
      RETURN
      END
