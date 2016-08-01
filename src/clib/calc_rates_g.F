#include "phys_const.def"
c=======================================================================
c//////////////////////  SUBROUTINE CALC_RATES_G  \\\\\\\\\\\\\\\\\\\\\\\
c
      subroutine calc_rates_g(
     &             ispecies,
     &             nratec, aye, temstart, temend,
     &             casebrates, threebody,
     &             uxyz, uaye, urho, utim,
     &             ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &             ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &             reHeII2a, reHeIIIa, brema, compa, 
     &             gammahacgs, gammaha,
     &             hyd01ka, h2k01a, vibha, rotha, rotla, 
     &             gpldl, gphdl, hdltea, hdlowa, ciecoa,
     &             gaHIa, gaH2a, gaHea, gaHpa, gaela, h2ltea, gasgra,
     &             k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &             k11a, k12a, k13a, k13dda, k14a, k15a, k16a, k17a,
     &             k18a, k19a, k20a, k21a, k22a, k23a,
     &             k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &             k57a, k58a, 
     &             ndratec, dtemstart, dtemend, h2dusta, 
     &             ncrna, ncrd1a, ncrd2a, ioutput
     &                     )
c
c  COMPUTE MULTISPECIES RATE LOOKUP TABLE
c
c  written by: Yu Zhang, Peter Anninos and Tom Abel
c  date:       
c  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
c  modified2: October, 1996 by GB; moved to AMR
c  modified3: February, 2000 by GB; added Toms new collisional rates
c  modified4: July, 2010 by Dan Reynolds; added case-B recombination rates
c  modified5: November, 2015 by Simon Glover; added several new processes
c
c  PURPOSE:
c    Construct tables for rate coefficients and some cooling functions 
c    versus the logarithm of temperature.
c
c  UNITS:
c    Most rate coefficients are computed in cgs units and then converted
c    into more convenient units.  This would usually be code units
c    (see units.src) but in many cases the natural code units are not
c    near unity (e.g. particles/box volume), so we include units and
c    constants from the equations themselves.  This is explained for
c    each of two types used: rate coefficients, cooling coefficients.
c    Note also that the densities in the rate equations are true, not
c    comoving densities and so we must include a factor of a^3.
c
c  PARAMETERS:
c    nnu      = number of frequency bins
c    everg    = ergs per eV (Rydberg constant) in cgs
c    evhz     = everg / h
c    mh       = mass of hydrogen atom
c    temstart = table start temperature (K)
c    temend   = table end temperature (K)
c    dtemstart = dust table start temperature (K)
c    dtemend  = dust table end temperature (K)
c    tevk     = conversion constant [eV] = [K] / tevk
c    ireco    = Recombination cooling flag (1 = on)
c    ibrco    = Brehmmstrahlung cooling flag (1 = on)
c    icico    = Collisional ionization flag (1 = on)
c    iceco    = Collisional excitation flag (1 = on)
c    iphoto_mol = Molecular photo-dissociation flag (1 = on)
c
c  INPUTS
c    nratec   = number of temperature bins in rate equations
c    ndratec  = number of dust temperature bins for H2+dust rate
c
c    tbase1   = code time units
c    xbase1   = code length units
c    kunit    = conversion factor xbase1**3/tbase1
c
c    casebrates = use case-B recombination rates for k2, k4 and k6
c
c  RATES:
c
c  Reaction rate numbering agrees with Abel et al. (1997).
c
c ------- 1:    HI    + e   -> HII   + 2e
c ------- 2:    HII   + e   -> H     + p
c ------- 3:    HeI   + e   -> HeII  + 2e
c ------- 4:    HeII  + e   -> HeI   + p
c ------- 5:    HeII  + e   -> HeIII + 2e
c ------- 6:    HeIII + e   -> HeII  + p
c ------- 7:    HI    + e   -> HM    + p
c ------- 8:    HM    + HI  -> H2I*  + e
c ------- 9:    HI    + HII -> H2II  + p
c ------- 10:   H2II  + HI  -> H2I*  + HII
c ------- 11:   H2I   + HII -> H2II  + H
c ------- 12:   H2I   + e   -> 2HI   + e
c ------- 13:   H2I   + H   -> 3H
c ------- 14:   HM    + e   -> HI    + 2e
c ------- 15:   HM    + HI  -> 2H    + e
c ------- 16:   HM    + HII -> 2HI
c ------- 17:   HM    + HII -> H2II  + e
c ------- 18:   H2II  + e   -> 2HI
c ------- 19:   H2II  + HM  -> HI    + H2I
c ------- 20:   unused
c ------- 21:   2H    + H2  -> H2I   + H2I
c ------- 22:   2H    + H   -> H2I   + H
c ------- 24:   HI    + p   -> HII   + e
c ------- 25:   HeII  + p   -> HeIII + e
c ------- 26:   HeI   + p   -> HeII  + e
c ------- 27:   HM    + p   -> HI    + e
c ------- 28:   H2II  + p   -> HI    + HII
c ------- 29:   H2I   + p   -> H2II  + e
c ------- 30:   H2II  + p   -> 2HII  + e
c ------- 31:   H2I   + p   -> 2HI
c ------- 50:   H+    + D   -> H     + D+
c ------- 51:   H     + D+  -> H+    + D
c ------- 52:   H2    + D+  -> HD    + H+
c ------- 53:   HD    + H+  -> H2    + D+
c ------- 54:   H2    + D   -> HD    + H
c ------- 55:   HD    + H   -> H2    + D
c ------- 56:   D     + H-  -> HD    + e-
c ------- 57:   D-    + H   -> HD    + e-  (included by multiplying 56 by 2)
c ------- 57:   HI    + HI  -> HII   + HI    + e
c ------- 58:   HI    + HeI -> HII   + HeI   + e
c
c ------- h2dust: 2H  + grain -> H2  + grain (eq 3.8 from Hollenbach & McKee 1989)
c
c ncrn, ncrd1, ncrd1: numerator and denominator terms in Eq 23 of Omukai ea 2000
c               for H2 formation heating.
c
c-----------------------------------------------------------------------
c
c
      implicit NONE
#include "grackle_fortran_types.def"
c
c  Arguments
c
      integer ispecies, nratec, casebrates, threebody, ndratec
      real*8  aye, temstart, temend, dtemstart, dtemend
      real*8  uxyz, uaye, urho, utim
      real*8  ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec),
     &        hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldl(nratec),
     &        gphdl(nratec), hdltea(nratec), hdlowa(nratec),
     &        gaHIa(nratec), gaH2a(nratec), gaHea(nratec), 
     &        gaHpa(nratec), gaela(nratec), h2ltea(nratec),
     &        gasgra(nratec), ciecoa(nratec),
     &        compa, gammaha, gammahacgs,
     &        k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &        k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &        k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &        k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &        k17a(nratec), k18a(nratec), k19a(nratec), k20a(nratec), 
     &        k21a(nratec), k23a(nratec),
     &        k22a(nratec), k50a(nratec), k51a(nratec), k52a(nratec), 
     &        k53a(nratec), k54a(nratec), k55a(nratec), k56a(nratec),
     &        k57a(nratec), k58a(nratec),
     &        k13dda(nratec, 14),
     &        h2dusta(nratec, ndratec),
     &        ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec)
      integer ioutput
c
c  Parameters
c
      integer ireco, ibrco, icico, iceco, iphoto_mol
      parameter (ireco = 1, ibrco = 1, icico = 1, iceco = 1, 
     &           iphoto_mol = 0)

      real*8 everg, evhz, tevk, mh, pi, tiny_cgs, dhuge, kb
      parameter (everg = ev2erg, evhz  = 2.41838d+14, 
     &           tevk = 1.1605d+4, mh = mass_h, 
     &           pi=pi_val, tiny_cgs = 1.0d-37, 
     &           dhuge = 1.0d+30, kb = kboltz)

c
c         Set various cutoff values in eV.
c
      real*8 e24,e25,e26,e27,e28a,e28b,e28c,e29a,e29b,e29c,
     &                 e30a,e30b
        PARAMETER(
     &    e24  = 13.6d0
     &   ,e25  = 54.4d0
     &   ,e26  = 24.6d0
     &   ,e27  = 0.755d0
     &   ,e28a = 2.65d0
     &   ,e28b = 11.27d0
     &   ,e28c = 21.0d0
     &   ,e29a = 15.42d0
     &   ,e29b = 16.5d0
     &   ,e29c = 17.7d0
     &   ,e30a = 30.0d0
     &   ,e30b = 70.0d0)

c  Locals
c
      integer i,j,idt
      real*8 logttt, ttt, tev, logtev, 
     &     xx, dum, tbase1, xbase1, kunit, coolunit,
     &     dbase1, dlogtem, kunit_3bdy, cierate,
     &     grain_coef, d_ttt, d_logttt, d_dlogtem,
     &     ttt2, d_ttt2, tk9, lambdaHI, lambdaHeII,
     &     lambdaHeIII
      real*8 tm, HDLR, HDLV, lt, t3, lt3
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////
c=======================================================================
c
c  Get conversion units
c
c     t/x/dbase1 is the number (z dependant) that converts from the
c       dimensionless code units to physical units.  Also, in the
c       code aye = 1 at z=zinit, so to convert the usual a (=1 at z=0)
c       to a~ (written in the code as aye), we use a = a~*[a] 
c
      tbase1 = utim
      xbase1 = uxyz/(aye*uaye)    ! uxyz is [x]*a     
      dbase1 = urho*(aye*uaye)**3 ! urho is [dens]/a^3
c
c  1) Set the dimensions of the (non-radiative) rate coefficients.  
c    Note that we have included the units that convert density to 
c    number density, so the rate equations should look like 
c    (in dimensionless units, hence the primes):
c
c       d(d0~)/dt~ = k~ * d1~ * d2~ / a~^3
c
c    where k~ is the dimenionless rate coefficients and d0-2~ are three
c     dimensionless densities (i.e. d = [dens]*d~) and a~ is the 
c     dimensionless expansion coefficient (see above).
c
c    rate eqn        : delta(n0)  = k  * n1        * n2        * dt     / a^3
c    rate eqn units  : [dens]/mh  = k  * [dens]/mh * [dens]/mh * [time] / [a]^3
c    rate eqn dimless: delta(n0~) = k~ * n1~       * n2~       * dt~    / a~^3
c    so: k = [k] * k~  where [k] = ( [a]^3 * mh ) / ( [dens] * [time] )  (~)
c    reminder: the number densities here are normalized with [dens] which
c              is not a constant (it has a factor a^3), so the number
c              densities must be converted from comoving to proper.
c
      kunit   = (uaye**3 * mh) / (dbase1 * tbase1)
      kunit_3bdy  = kunit * (uaye**3 * mh) / dbase1
c
c  2) Set the dimension of the cooling coefficients (including constants)
c     (this equation has a rho because e is the specific energy, not
c      energy/unit volume).
c       delta(e)  = L     * n1        * n2        * dt     / dens   / a^3
c       [e]       = L     * [dens]/mh * [dens]/mh * [time] / [dens] / [a]^3
c       delta(e~) = L~    * n1~       * n2~       * dt~    / dens~  / a~^3 [~]
c     so L = [L] * L~ where [L] = [e] * mh**2 * [a]^3 / ([dens] * [time]) [~]
c       but [e] = ([a]*[x])**2 / [time]**2 and ([a] = 1 / (1 + zri) )
c      [L] = ([a]**5 * [x]**2 * mh**2) / ([dens] * [time]**3)
c
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
c
c    Note: some of the coffiecients have only one power of n.  These
c          do not have the /a^3 factor, also they have units
c          [L1] = ([a]**2 * [x]**2 * mh) / [time]**3
c               = [L] * [dens] * [a]**3 / mh
c          This is done through the dom variable in cool.src
c         (some have three powers of n and they are different by the
c          reciprocal of the above factor multiplying [L]).
c
c  3) the units for the radiative rate coefficients is just 1/[time]
c
c  4) Energy transfer from gas to dust grains, following equation 2.15
c     of Hollenbach & McKee (1989).
      grain_coef = 1.2d-31 * 1.0d3**(-0.5d0)
c
c  Compute log spacing in temperature
c
      ttt    = temstart
      logttt = log(ttt)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1, 8)
c
c  Compute log spacing in dust temperature
c
      d_ttt     = dtemstart
      d_logttt  = log(d_ttt)
      d_dlogtem = (log(dtemend) - log(dtemstart))/real(ndratec-1, 8)

      if (ispecies .gt. 0) then

c
c  Initialize constants to tiny
c
      do i = 1, nratec
c
        k1a (i) = tiny
        k2a (i) = tiny
        k3a (i) = tiny
        k4a (i) = tiny
        k5a (i) = tiny
        k6a (i) = tiny
        k7a (i) = tiny
        k8a (i) = tiny
        k9a (i) = tiny
        k10a(i) = tiny
        k11a(i) = tiny
        k12a(i) = tiny
        k13a(i) = tiny
        k14a(i) = tiny
        k15a(i) = tiny
        k16a(i) = tiny
        k17a(i) = tiny
        k18a(i) = tiny
        k19a(i) = tiny
        k20a(i) = tiny
        k21a(i) = tiny
        k22a(i) = tiny
        k50a(i) = tiny
        k51a(i) = tiny
        k52a(i) = tiny
        k53a(i) = tiny
        k54a(i) = tiny
        k55a(i) = tiny
        k56a(i) = tiny
        k57a(i) = tiny
        k58a(i) = tiny
c
        do j = 1, 14
           k13dda(i,j) = tiny
        enddo
c
        do j = 1, ndratec
           h2dusta(i,j) = tiny
        enddo
c
        ncrna(i)  = tiny
        ncrd1a(i) = tiny
        ncrd2a(i) = tiny
c
        ceHIa(i)    = tiny
        ceHeIa(i)   = tiny
        ceHeIIa(i)  = tiny
        ciHIa(i)    = tiny
        ciHeIa(i)   = tiny
        ciHeISa(i)  = tiny
        ciHeIIa(i)  = tiny
        reHIIa(i)   = tiny
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
        brema(i)    = tiny
        compa       = tiny
c
        hyd01ka(i)  = tiny
        h2k01a(i)   = tiny
        vibha(i)    = tiny
        rotha(i)    = tiny
        rotla(i)    = tiny
        hdltea(i)   = tiny
        hdlowa(i)   = tiny
        ciecoa(i)   = tiny
c
        gaHIa(i)    = tiny
        gaH2a(i)    = tiny
        gaHea(i)    = tiny
        gaHpa(i)    = tiny
        gaela(i)    = tiny
        h2ltea(i)   = tiny
        gasgra(i)   = tiny
      enddo
c
c  Fill in tables over the range temstart to temend
c
c -------------------------------------------------
c  1) rate coefficients (excluding external radiation field)
c
      do i = 1, nratec
c
c       Compute temperature of this bin (in eV)
c
        logttt = log(temstart) + real(i-1, 8)*dlogtem
        ttt = exp(logttt)
        ttt2 = ttt / 1.0d2
        tev = ttt/tevk
        logtev = log(tev)
c
c
c       Call Tom Abels routine (from his web page, Feb/2000)
c
        call coll_rates_g(ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                  k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                  k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                  k16a(i), k17a(i), k18a(i), k19a(i), k23a(i),
     &                  kunit, casebrates)
c
c       Compute the coefficients needed for the density-dependent collision 
c       H2 dissociation rate. Our source for this (Martin et al, 1996, ApJ, 461, 265) 
c       distinguishes between two processes, direct collisional dissociation (idt=0) and
c       dissociative tunnelling (idt=1). We need to account for both processes, so call 
c       the routine twice, once for each process.
c
c       For each process, the routine returns seven coefficients that are later used
c       to compute the log10 of the corresponding rate.
c
        idt = 0
        call colh2diss_g(ttt, idt, k13dda(i,1), k13dda(i,2), 
     &                 k13dda(i,3), k13dda(i,4), k13dda(i,5), 
     &                 k13dda(i,6), k13dda(i,7))
c       Normalize by dividing by kunit
        k13dda(i,1) = k13dda(i,1) - log10(kunit)
c
        idt = 1
        call colh2diss_g(ttt, idt, k13dda(i,8), k13dda(i,9), 
     &                 k13dda(i,10), k13dda(i,11), k13dda(i,12), 
     &                 k13dda(i,13), k13dda(i,14))
c       Normalize by dividing by kunit
        k13dda(i,8) = k13dda(i,8) - log10(kunit)
c
c       ------ 3-body H2 rate ----
c       The first bit is my fit to A.E. Orel 1987, J.Chem.Phys., 87, 
c       314. I then match it to the T^-1 of Palla etal (1983) 
c       Which is then 4 times smaller than the Palla rate ! Thats
c       molecule uncertainties ! :-)
c
c Our varying threebody rates from Simon, with corresponding 
c collisional dissociation rates
        if(threebody.eq.0)then
          if (ttt .le. 300.d0) then
             k22a(i) = 1.3d-32 * (ttt/300.d0)**(-0.38d0)
          else
             k22a(i) = 1.3d-32 * (ttt/300.d0)**(-1.d0)
          endif
        elseif(threebody.eq.1)then
c PSS83 three-body rate
          k22a(i) = 5.5d-29 / ttt
c Inverse of PSS83 three-body rate
          k13a(i) = (5.24d-7 / ttt**0.485d0) 
     &                   * exp(-5.2d4 / ttt)
        elseif(threebody.eq.2)then
c CW83 three-body rate
          k22a(i) = 8.8d-33
c Inverse of CW83 rate
          k13a(i) = 8.4d-11 * ttt**0.515d0 * exp(-5.2d4 / ttt)
        elseif(threebody.eq.3)then
c FH07 three-body rate
          k22a(i) = 1.44d-26 / ttt**1.54d0
c Rate from JGC67, used by FH07 to construct their three-body rate
          k13a(i) = (1.38d-4 / ttt**1.025d0)
     &                   * exp(-5.2d4 / ttt)
        elseif(threebody.eq.4)then
c Rate from Glover (2008), derived by detailed balance from MSM96
          k22a(i) = 7.7d-31 / ttt**0.464d0
c High density limiting rate from MSM96
          k13a(i) = 1d1**(-178.4239d0 - 68.42243d0 * log10(ttt) 
     &                   + 43.20243d0 * log10(ttt)**2
     &                   - 4.633167d0 * log10(ttt)**3 
     &                   + 69.70086d0 * log10(1d0 + 40870.38d0 / ttt)
     &                   - (23705.7d0 / ttt))
        elseif(threebody.eq.5)then
c Rate from Forrey (2013)
          k22a(i) = (6d-32 / ttt**0.25d0) + (2d-31 / ttt**0.5d0)
c Inferred from detailed balance from Forrey rate
          if (ttt .le. 3000.d0) then
            k13a(i) = 2.4d-8 * exp(-5.2d4/ttt)
          else
            k13a(i) = 2.2d-6 * ttt**(-0.565) * exp(-5.2d4/ttt)
          endif
        else
            write(0,*) "I didn't understand your threebody rate."
            stop
        endif
        k13a(i) = k13a(i) / kunit
        k22a(i) = k22a(i) / kunit_3bdy
        
        k21a(i) = 2.8d-31 * (ttt**(-0.6d0)) / kunit_3bdy
c        write(6,*),'calc_rates_g: check: ', ttt, k2a(i), k22a(i)
c
c       ------ Deuterium rates -----
c
c       50) H+  +  D  -> H  +  D+
c       51) H   +  D+ -> H+ +  D
c       52) H2  +  D+ -> HD +  H+
c       53) HD  +  H+ -> H2 +  D+
c       54) H2  +  D  -> HD +  H
c       55) HD  +  H  -> H2 +  D
c       56) D   +  H- -> HD +  e-
c      [57) D-  +  H  -> HD +  e-]  included by multiply 56 by 2
c
        if (ttt .le. 2d5) then
c Fit from Savin (2002), valid for T < 2e5 K; at higher T, we extrapolate
          k50a(i) = (2.0d-10 * ttt**0.402d0 * exp(-3.71d1/ttt)
     &             - 3.31d-17 * ttt**1.48d0) / kunit
        else
          k50a(i) = 2.5d-8 * (ttt/2.0d5)**0.402d0 / kunit
        endif
c Fit from Savin (2002), valid for T < 2e5 K
        k51a(i) = (2.06d-10 * ttt**0.396d0 * exp(-3.30d1/ttt)
     &            + 2.03d-9 * ttt**(-0.332d0)) / kunit
c k52, k53: fits by Galli & Palla (2002) to calculations by Gerlich (1982)
c For k52, we use a fixed value at T > 1e4 K to avoid numerical problems with
c the fitting function. Note that at these temperatures, we do not expect the
c reaction to be important
        if (ttt. le. 1d4) then
          k52a(i) = 1.0d-9 * (0.417 + 0.846 * log10(ttt) 
     &            - 0.137 * log10(ttt)**2)   / kunit
        else
          k52a(i) = 1.609d-9 / kunit
        endif
        k53a(i) = 1.1d-9 *exp(-4.88d2/ttt) / kunit
c Fit by Clark et al (2011), based on data from in Mielke et al (2003)
        if (ttt .le. 2d3) then
          k54a(i) = 1d1**(-5.64737d1 + 5.88886 * log10(ttt)
     $            + 7.19692d0  * log10(ttt)**2
     $            + 2.25069d0  * log10(ttt)**3
     $            - 2.16903d0  * log10(ttt)**4
     $            + 3.17887d-1 * log10(ttt)**5)
        else
          k54a(i) = 3.17d-10 * exp(-5.207d3 / ttt)
        endif
c From Galli & Palla (2002), based on Shavitt (1959). Fit modified at low
c temperature, as recommended by Ripamonti (2007) and McGreer & Bryan (2008)
c to avoid producing an anomalously large rate coefficient
        if (ttt .le. 2d2) then
          k55a(i) = 1.08d-22 / kunit
        else
          k55a(i) = 5.25d-11 * exp(-4.43d3/ttt + 1.739d5/ttt**2) / kunit
        endif
c Same as for H- + H; measurements by Miller et al (2012) suggest
c that there is no significant isotope effect for this reaction
        k56a(i) = k8a(i)
c
c       ------ New H ionization rates -----
c
c       57) H  +  H  ->  H+ + H  + e
c       58) H  +  He ->  H+ + He + e
c
c       The rate coefficients for reactions 57 and 58 are from 
c       Lenzuni, Chernoff & Salpeter (1991). The value for reaction 57 is based
c       on the experimental cross-sections of Gealy & van Zyl (1987); that for
c       reaction 58 is based on cross-sections from van Zyl, Le & Amme (1981).
c
        if (ttt .gt. 3d3) then
          k57a(i) = 1.2d-17  * ttt**1.2d0 * dexp(-1.578d5 / ttt) / kunit
          k58a(i) = 1.75d-17 * ttt**1.3d0 * dexp(-1.578d5 / ttt) / kunit
        else
          k57a(i) = tiny
          k58a(i) = tiny
        endif
c
c     H2 formation on dust grains.
c     Loop over dust temperature.
c
        do j = 1, ndratec
c
c     Compute dust temperature of this bin
c
           d_logttt = log(dtemstart) + real(j-1, 8)*d_dlogtem
           d_ttt = exp(d_logttt)
           d_ttt2 = d_ttt / 1.0d2
c
c     k23 from Omukai (2000).
c
           h2dusta(i,j) = 6.0d-17 * (ttt / 300.0d0)**(0.5d0) *
     &          ((1.0d0 + exp(7.5d2 * ((1.0d0 / 75.0d0) - 
     &                                (1.0d0 / d_ttt))))**(-1.0d0)) *
     &          ((1.0d0 + (4.0d-2 * (ttt + d_ttt)**(0.5d0)) + 
     &           (2.0d-3 * ttt) + 
     &           (8.0d-6 * ttt**(2.0d0)))**(-1.0d0)) / kunit
c
c     Equation 3.8 from Hollenbach & McKee (1979).
c
c           h2dusta(i,j) = 3.0d-17 * ttt2**0.5d0 /
c     &          (1.0d0 + 
c     &           (0.4d0 * (ttt2 + d_ttt2)**0.5d0) +
c     &           (0.2d0 * ttt2) +
c     &           (8.0d-2 * ttt2**2.0d0)) / kunit
c
        enddo
c
c       H2 formation heating terms from Equation 23 of Omukai (2000).
c
        ncrna(i)  = 1.0d6 * (ttt**(-0.5d0))
        ncrd1a(i) = 1.6d0 * exp(-(400.0d0 / ttt)**(2.d0))
        ncrd2a(i) = 1.4d0 * exp(-12000.0d0 / (ttt + 1200.0d0))
c
      enddo
c
c  Write out cgs rate coefficients
c
c      open(10, file='k1-6.out', status='unknown')
c      do i = 1, nratec
c         write(10,1010) exp(log(temstart) + real(i-1,8)*dlogtem),
c     &               k1a(i),k2a(i),k3a(i),k4a(i),k5a(i),k6a(i)
c      enddo
c 1010 format(7(1pe11.4))
c      close(10)
c
c -------------------------------------------------
c  2) Cooling/heating rates (excluding external radiation field)
c
      do i = 1, nratec
c
        logttt = log(temstart) + real(i-1, 8)*dlogtem
        ttt = exp(logttt)
c
c    a) Collisional excitations (Black 1981; Cen 1992)
c
      if ( iceco .eq. 1 ) then
        ceHIa(i)    = 7.5d-19*exp(-min(log(dhuge),118348.d0/ttt))
     &              /(1.d0+sqrt(ttt/1.0d5)) / coolunit
        ceHeIa(i)   = 9.1d-27*exp(-min(log(dhuge),13179.d0/ttt))
     &              *ttt**(-0.1687d0)/(1.d0+sqrt(ttt/1.0d5)) / coolunit
        ceHeIIa(i)  = 5.54d-17*exp(-min(log(dhuge),473638.d0/ttt))
     &              *ttt**(-0.397d0)/(1.d0+sqrt(ttt/1.0d5)) / coolunit
      else
        ceHIa(i)    = tiny
        ceHeIa(i)   = tiny
        ceHeIIa(i)  = tiny
      endif
c
c    b) Collisional ionizations (Cen 1992 or Abel 1996)
c
      if ( icico .eq. 1 ) then
c
c        ciHIa(i)    = 1.27d-21*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),157809.1/ttt)) / coolunit
c        ciHeIa(i)   = 9.38d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),285335.4/ttt)) / coolunit
c        ciHeIIa(i)  = 4.95d-22*sqrt(ttt)/(1.+sqrt(ttt/1.0d5))
c     &              *exp(-min(log(dhuge),631515.0/ttt)) / coolunit
c
        ciHeISa(i)  = 5.01d-27*(ttt)**(-0.1687d0)/(1.d0+sqrt(ttt/1.0d5))
     &              *exp(-min(log(dhuge),55338.d0/ttt)) / coolunit
c
c       Collisional ionizations (polynomial fits from Tom Abel)
c
        ciHIa(i)    = 2.18d-11*k1a(i)*kunit / coolunit
        ciHeIa(i)   = 3.94d-11*k3a(i)*kunit / coolunit
        ciHeIIa(i)  = 8.72d-11*k5a(i)*kunit / coolunit
      else
        ciHIa(i)    = tiny
        ciHeIa(i)   = tiny
        ciHeIIa(i)  = tiny
        ciHeISa(i)  = tiny
      endif
c
c    c) Recombinations (from Hui & Gnedin 1997 except where noted)
c
      if ( ireco .eq. 1 ) then
        lambdaHI    = 2.d0 * 157807.d0 / ttt
        lambdaHeII  = 2.d0 * 285335.d0 / ttt
        lambdaHeIII = 2.d0 * 631515.d0 / ttt
        if (casebrates .eq. 1) then
          reHIIa(i)   = 3.435d-30 * ttt * lambdaHI**1.970
     &                / (1.d0 + (lambdaHI/2.25d0)**0.376)**3.720 
     &                / coolunit
          reHeII1a(i)  = 1.26d-14 * kb * ttt * lambdaHeII**0.75
     &                 / coolunit
          reHeIIIa(i) = 8.d0 * 3.435d-30 * ttt * lambdaHeIII**1.970
     &                / (1.d0 + (lambdaHeIII/2.25d0)**0.376)**3.720 
     &                / coolunit
        else
          reHIIa(i)   = 1.778d-29 * ttt * lambdaHI**1.965
     &                / (1.d0 + (lambdaHI/0.541d0)**0.502)**2.697 
     &                / coolunit
          reHeII1a(i) = 3d-14 * kb * ttt * lambdaHeII**0.654d0 
     &                / coolunit
          reHeIIIa(i) = 8.d0 * 1.778d-29 * ttt * lambdaHeIII**1.965
     &                / (1.d0 + (lambdaHeIII/0.541d0)**0.502)**2.697 
     &                / coolunit
        endif
c Dielectronic recombination: rate from Cen (1992). Same rate regardless
c of whether we have case A or case B.
        reHeII2a(i) = 1.24d-13*ttt**(-1.5d0)
     &              *exp(-min(log(dhuge),470000.d0/ttt))
     &              *(1.d0+0.3d0*exp(-min(log(dhuge),94000.d0/ttt))) 
     &              / coolunit
      else
        reHIIa(i)   = tiny
        reHeII1a(i) = tiny
        reHeII2a(i) = tiny
        reHeIIIa(i) = tiny
      endif
c
c    d) Bremsstrahlung (Black 1981)(Spitzer & Hart 1979)
c
      if ( ibrco .eq. 1 ) then
        brema(i)    = 1.43d-27*sqrt(ttt)
     &              *(1.1d0+0.34d0*exp(-(5.5d0-log10(ttt))**2/3.d0)) / 
     &               coolunit
      else
        brema(i)    = tiny
      endif
c
c         Bremsstrahlung (Shapiro & Kang 1987)(Spitzer 1978)
c
c        if(ttt .lt. 1.0*3.2d5) gaunt = 0.79464 + 0.1243*log10(ttt/1.)
c        if(ttt .ge. 1.0*3.2d5) gaunt = 2.13164 - 0.1240*log10(ttt/1.)
c        brem1a(i) = 1.426d-27*sqrt(ttt)*gaunt / coolunit
c        if(ttt .lt. 4.0*3.2d5) gaunt = 0.79464 + 0.1243*log10(ttt/4.)
c        if(ttt .ge. 4.0*3.2d5) gaunt = 2.13164 - 0.1240*log10(ttt/4.)
c        brem2a(i) = 1.426d-27*sqrt(ttt)*gaunt / coolunit
c
c    e) Molecular hydrogen cooling
c
c       (Which one is used is controlled by a flag in cool1d_mulit_g.F)
c
c    e part1) - Lepp & Shull rates
c
        xx = dlog10(ttt/1.0d4)
        vibha  (i) = 1.1d-18*exp(-min(log(dhuge),6744.d0/ttt)) 
     &               / coolunit
c
        if( ttt .gt. 1635.d0) then
           dum = 1.0d-12*sqrt(ttt)*exp(-1000.d0/ttt)
        else
           dum = 1.4d-13*exp((ttt/125.d0) - (ttt/577.d0)**2)
        endif
        hyd01ka(i) = dum*exp(-min(log(dhuge),
     &                            8.152d-13/(kboltz*ttt) )) 
     &               / coolunit
c
        dum = 8.152d-13*(4.2d0/(kboltz*(ttt+1190.d0)) +
     &       1.d0/(kboltz*ttt))
        h2k01a(i) = 1.45d-12*sqrt(ttt)*exp(-min(log(dhuge),dum)) 
     &              / coolunit
c
        if(ttt .gt. 4031.d0)then
          rotla(i) = 1.38d-22*exp(-9243.d0/ttt) / coolunit
        else
          rotla(i) = 10.d0**(-22.9d0 - 0.553d0*xx - 1.148d0*xx*xx) / 
     &               coolunit
        endif
c
        if(ttt .gt. 1087.d0)then
          rotha(i) = 3.9d-19*exp(-6118.d0/ttt) / coolunit
        else
          rotha(i) = 10.d0**(-19.24d0 + 0.474d0*xx - 1.247d0*xx*xx) / 
     &               coolunit
        endif
c
c     e part2) - Galli and Palla (1999) rates as fit by Tom Abel
c
        tm = max(ttt, 13.0d0)       ! no cooling below 13 Kelvin ...
        tm = min(tm, 1.d5)      ! fixes numerics
        lt = dlog10(tm)
c     low density limit from Galli and Palla
        gpldl(i)  = 10.d0**(-103.d0+97.59d0*lt-48.05d0*lt**2 
     &       +10.8d0*lt*lt*lt-0.9032d0*lt*lt*lt*lt) / coolunit
c     high density limit from HM79 (typo corrected Aug 30/2007)
        t3 = tm/1000.d0
        HDLR = ((9.5e-22*t3**3.76d0)/(1.d0+0.12d0*t3**2.1d0)*
     &            exp(-(0.13d0/t3)**3)+3.d-24*exp(-0.51d0/t3))
        HDLV = (6.7d-19*exp(-5.86d0/t3) + 1.6d-18*exp(-11.7d0/t3))
        gphdl(i)  = (HDLR + HDLV) / coolunit
c
c     e part 3) - Low density rates from Glover & Abel (2008)
c
c  Excitation by HI
c
      tm  = max(ttt, 10.0d0)
      tm  = min(tm, 1.d4)
      lt3 = log10(tm / 1.d3)  

      if (tm .lt. 1e2) then
        gaHIa(i) = 10.d0**(-16.818342d0
     &         + 37.383713d0 * lt3
     &         + 58.145166d0 * lt3**2
     &         + 48.656103d0 * lt3**3
     &         + 20.159831d0 * lt3**4
     &         + 3.8479610d0 * lt3**5) / coolunit
      elseif (tm .lt. 1e3) then
        gaHIa(i) = 10.d0**(-24.311209d0
     &         + 3.5692468d0 * lt3
     &         - 11.332860d0 * lt3**2
     &         - 27.850082d0 * lt3**3
     &         - 21.328264d0 * lt3**4
     &         - 4.2519023d0 * lt3**5) / coolunit
      else
        gaHIa(i) = 10.d0**(-24.311209d0
     &         + 4.6450521d0 * lt3
     &         - 3.7209846d0 * lt3**2
     &         + 5.9369081d0 * lt3**3
     &         - 5.5108047d0 * lt3**4
     &         + 1.5538288d0 * lt3**5) / coolunit
      endif
c
c Excitation by H2
c
      gaH2a(i) = 10.d0**(-23.962112d0
     &        + 2.09433740d0  * lt3
     &        - 0.77151436d0  * lt3**2
     &        + 0.43693353d0  * lt3**3
     &        - 0.14913216d0  * lt3**4
     &        - 0.033638326d0 * lt3**5) / coolunit
c
c Excitation by He
c
      gaHea(i) = 10.d0**(-23.689237d0
     &        + 2.1892372d0  * lt3
     &        - 0.81520438d0 * lt3**2
     &        + 0.29036281d0 * lt3**3
     &        - 0.16596184d0 * lt3**4
     &        + 0.19191375d0 * lt3**5) / coolunit
c
c Excitation by H+
c
c Revised rate, using new collisional excitation rates from
c Honvault et al (2011, Phys. Rev. Lett., 107, 023201) and
c Honvault et al (2012, Phys. Rev. Lett., 108, 109903).
c
      gaHpa(i) = 10.d0**(-22.089523
     &         + 1.5714711   * lt3
     &         + 0.015391166 * lt3**2
     &         - 0.23619985  * lt3**3
     &         - 0.51002221  * lt3**4
     &         + 0.32168730  * lt3**5) / coolunit
c
c Excitation by electrons
c
c Revised rate, based on data from  Yoon et al (2008,
c J. Phys. Chem. Ref. Data, 37, 913).
c
      if (tm .lt. 100) then
        gaela(i) = 0d0
      elseif (tm .lt. 500) then
        gaela(i) = 10.d0**(-21.928796d0
     &           + 16.815730d0 * lt3
     &           + 96.743155d0 * lt3**2
     &           + 343.19180d0 * lt3**3
     &           + 734.71651d0 * lt3**4
     &           + 983.67576d0 * lt3**5
     &           + 801.81247d0 * lt3**6
     &           + 364.14446d0 * lt3**7
     &           + 70.609154d0 * lt3**8) / coolunit
      else
        gaela(i) = 10.d0**(-22.921189d0
     &           + 1.6802758d0  * lt3
     &           + 0.93310622d0 * lt3**2
     &           + 4.0406627d0  * lt3**3
     &           - 4.7274036d0  * lt3**4
     &           - 8.8077017d0  * lt3**5
     &           + 8.9167183d0  * lt3**6
     &           + 6.4380698d0  * lt3**7
     &           - 6.3701156d0  * lt3**8) / coolunit
      endif
c     e part 4) - New fit to LTE rate from Glover (2015, MNRAS, 451, 2082)
      if (tm .lt. 1d2) then
c Crude extrapolation, but don't expect H2 cooling to be significant
c at these temperatures
        h2ltea(i) = 7.d-27 * tm**1.5 * exp(-512./tm) / coolunit
      else
        h2ltea(i) = 10.d0**(-20.584225d0
     &            + 5.0194035d0 * lt3
     &            - 1.5738805d0 * lt3**2
     &            - 4.7155769d0 * lt3**3
     &            + 2.4714161d0 * lt3**4
     &            + 5.4710750d0 * lt3**5
     &            - 3.9467356d0 * lt3**6
     &            - 2.2148338d0 * lt3**7
     &            + 1.8161874d0 * lt3**8) / coolunit
      endif
c
c     f) HD cooling
c
c        HD Cooling Function (ergs cm3 /s)
c
c Fit from Coppola et al 2011 LTE (ergs/s) -> hdlte (ergs cm3/s)
c
      tm  = max(ttt, 10.0d0)
      tm  = min(tm, 3.d4)
      hdltea(i) = -55.5725d0 + 56.649d0 * dlog10(tm) -
     &            37.9102d0  * (dlog10(tm))**2 + 
     &            12.698d0   * (dlog10(tm))**3 - 
     &            2.02424d0  * (dlog10(tm))**4 + 
     &            0.122393d0 * (dlog10(tm))**5
      hdltea(i) = (10.d0**min(hdltea(i), 0.0d0)) / coolunit
c
c Rate based on Wrathmall, Gusdorf & Flower (2007) collisional excitation rates
c for HD-H collisions
c
      tm = max(ttt, 1.d1)
      tm = min(tm, 6.d3)

      hdlowa(i) = -23.175780d0  + 1.5035261d0 * dlog10(tm/1.d3)
     &          + 0.40871403d0  * (dlog10(tm/1.d3))**2
     &          + 0.17849311d0  * (dlog10(tm/1.d3))**3
     &          - 0.077291388d0 * (dlog10(tm/1.d3))**4
     &          + 0.10031326d0  * (dlog10(tm/1.d3))**5
      hdlowa(i) = (10.d0**hdlowa(i)) / coolunit
c
c     g) CIE cooling as discussed in Ripamonti & Abel (2003)
c        cie_thin_cooling_rate_g returns rate as erg s^-1 cm^3 g^2 so
c         divide my mh^2 to get usual erg s^-1 cm^3, which is then
c         converted to code units with the conversion factor coolunit.
c         To get cooling rate
c     
      call cie_thin_cooling_rate_g(ttt, cierate)
c      ciecoa(i) = cierate/mh/mh/coolunit
c
c     Matt Turk: cie_thin_cooling actually returns
c     erg per second times cm^3 per gram per (gram in H2 molecules)
c     To reproduce eq 5 in RA04, divide c_t_c_r by mh, so to get
c     erg/s cm^3 we multiply by mh.  We can either divide by 2 for 
c     the H2 mass->number or in cool1d_multi_g.F, but we might
c     as well do it here.
c
      ciecoa(i) = cierate * (mh/2.0d0) / coolunit
c      write(6,*), "thin: ", ttt, cierate * 1e14*mh*2
c
c     h) Energy transfer from gas to dust grains, following 
c        equation 2.15 of Hollenbach & McKee (1989)
      gasgra(i) = grain_coef * (ttt**0.5d0) * 
     &     (1.0d0 - (0.8d0 * exp(-75.0d0 / ttt))) / coolunit
c
      enddo

      endif
c
c     i) Compton cooling (Peebles 1971)
c
      compa     = 5.65d-36 / coolunit
c
c     j) Photoelectric heating by UV-irradiated dust (Wolfire 1995)
c        Default is 8.5e-26 for epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3)
c
      gammaha   = gammahacgs / coolunit

c
c     Write out cooling rate data
c

      if (ioutput.eq.1) then
        open(10, file='cool_rates.out', status='unknown')

        do i=1, nratec
           logttt = log(temstart) + real(i-1, 8)*dlogtem
           ttt = exp(logttt)
           write(10,1000) ttt, ceHIa(i), ceHeIa(i), ceHeIIa(i),
     &                    ciHIa(i), ciHeIa(i), ciHeISa(i),
     &                    ciHeIIa(i), reHIIa(i), reHeII1a(i),
     &                    reHeII2a(i), reHeIIIa(i), brema(i),
     &                    compa       
 1000      format(1p,30(e10.3,1x))
        enddo

        close(10)

c
c     Write out species rate data
c

        open(10, file='rates.out', status='unknown')

        do i=1, nratec
           logttt = log(temstart) + real(i-1, 8)*dlogtem
           ttt = exp(logttt)
           write(10,1000) ttt, k1a(i), k2a(i), k3a(i), k4a(i), k5a(i),
     &                    k6a(i), k7a(i), k8a(i), k9a(i), k10a(i),
     &                    k11a(i), k12a(i), k13a(i), k14a(i), k15a(i), 
     &                    k16a(i), k17a(i), k18a(i), k19a(i), k22a(i)
        enddo

       close(10)

      endif

      return
      end
