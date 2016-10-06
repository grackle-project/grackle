#include "phys_const.def"

!=======================================================================
!///////////////////  SUBROUTINE SOLVE_RATE_COOL_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_rate_cool_g(icool, d, e, u, v, w, de,
     &                HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, iexpand, 
     &                ispecies, imetal, imcool, idust, idim,
     &                is, js, ks, ie, je, ke, ih2co, ipiht, igammah,
     &                dt, aye, temstart, temend, 
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh, dtoh, z_solar,
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a,
     &                k16a, k17a, k18a, k19a, k22a,
     &                k24, k25, k26, k27, k28, k29, k30, k31,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a,
     &                k57a, k58a,
     &                ndratec, dtemstart, dtemend, h2dusta,
     &                ncrna, ncrd1a, ncrd2a,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp, piHI, piHeI, piHeII,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla, 
     &                gpldla, gphdla, hdltea, hdlowa,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela, 
     &                h2ltea, gasgra, 
     &                iradshield, avgsighi, avgsighei, avgsigheii,
     &                iradtrans, iradcoupled, iradstep,
     &                irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                photogamma,
     &                ierr,
     &                ih2optical, iciecool, ithreebody, ciecoa, 
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat)

!
!  SOLVE MULTI-SPECIES RATE EQUATIONS AND RADIATIVE COOLING
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1:  January, 1996 by Greg Bryan; converted to KRONOS
!  modified2:  October, 1996 by GB; adapted to AMR
!  modified3:  May,     1999 by GB and Tom Abel, 3bodyH2, solver, HD
!  modified4:  June,    2005 by GB to solve rate & cool at same time
!  modified5:  April,   2009 by JHW to include radiative transfer
!  modified6:  September, 2009 by BDS to include cloudy cooling
!
!  PURPOSE:
!    Solve the multi-species rate and cool equations.
!
!  INPUTS:
!    icool    - flag to update energy from radiative cooling
!    in,jn,kn - dimensions of 3D fields
!
!    d        - total density field
!    de       - electron density field
!    HI,HII   - H density fields (neutral & ionized)
!    HeI/II/III - He density fields
!    DI/II    - D density fields (neutral & ionized)
!    HDI      - neutral HD molecule density field
!    HM       - H- density field
!    H2I      - H_2 (molecular H) density field
!    H2II     - H_2+ density field
!    kph*     - photoionization fields
!    gamma*   - photoheating fields
!
!    is,ie    - start and end indices of active region (zero based)
!    iexpand  - comoving coordinates flag (0 = off, 1 = on)
!    idim     - dimensionality (rank) of problem
!    ispecies - chemistry module (1 - H/He only, 2 - molecular H, 3 - D) 
!    imetal   - flag if metal field is active (0 = no, 1 = yes)
!    imcool   - flag if there is metal cooling
!    idust    - flag for H2 formation on dust grains
!    ih2co    - flag to include H2 cooling (1 = on, 0 = off)
!    ipiht    - flag to include photoionization heating (1 = on, 0 = off)
!    iradshield - flag for approximate self-shielding of UV background
!    avgsighi   - spectrum averaged ionization crs for HI for use with shielding
!    avgsighei  - spectrum averaged ionization crs for HeI for use with shielding
!    avgsigheii - spectrum averaged ionization crs for HeII for use with shielding
!    iradtrans - flag to include radiative transfer (1 = on, 0 = off)
!    iradcoupled - flag to indicate coupled radiative transfer
!    iradstep  - flag to indicate intermediate coupled radiative transfer timestep
!    irt_honly - flag to indicate applying RT ionization and heating to HI only

!    fh       - Hydrogen mass fraction (typically 0.76)
!    dtoh     - Deuterium to H mass ratio
!    z_solar  - Solar metal mass fraction
!    dt       - timestep to integrate over
!    aye      - expansion factor (in code units)
!
!    utim     - time units (i.e. code units to CGS conversion factor)
!    uaye     - expansion factor conversion factor (uaye = 1/(1+zinit))
!    urho     - density units
!    uxyz     - length units
!    utem     - temperature(-like) units
!
!    temstart, temend - start and end of temperature range for rate table
!    nratec   - dimensions of chemical rate arrays (functions of temperature)
!    dtemstart, dtemend - start and end of dust temperature range
!    ndratec  - extra dimension for H2 formation on dust rate (dust temperature)
!
!    icmbTfloor - flag to include temperature floor from cmb
!    iClHeat    - flag to include cloudy heating
!    priGridRank - rank of cloudy primordial cooling data grid
!    priGridDim  - array containing dimensions of cloudy primordial data
!    priPar1, priPar2, priPar3 - arrays containing primordial grid parameter values
!    priDataSize - total size of flattened 1D primordial cooling data array
!    priCooling  - primordial cooling data
!    priHeating  - primordial heating data
!    priMMW      - primordial mmw data
!    metGridRank - rank of cloudy metal cooling data grid
!    metGridDim  - array containing dimensions of cloudy metal data
!    metPar1, metPar2, metPar3 - arrays containing metal grid parameter values
!    metDataSize - total size of flattened 1D metal cooling data array
!    metCooling  - metal cooling data
!    metHeating  - metal heating data
!    iVheat      - flag for using volumetric heating rate
!    iMheat      - flag for using specific heating rate
!    Vheat       - array of volumetric heating rates
!    Mheat       - array of specific heating rates
!
!  OUTPUTS:
!    update chemical rate densities (HI, HII, etc)
!
!  PARAMETERS:
!    itmax   - maximum allowed sub-cycle iterations
!    mh      - H mass in cgs units
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"
#ifdef _OPENMP
#include "omp_lib.h"
#endif

!  General Arguments

      integer icool, in, jn, kn, is, js, ks, ie, je, ke, nratec, 
     &        iexpand, ih2co, ipiht, ispecies, imetal, idim,
     &        ierr, imcool, idust, 
     &        igammah, ih2optical, iciecool, ithreebody,
     &        ndratec, clnew, iVheat, iMheat, iradshield,
     &        iradtrans, iradcoupled, iradstep, irt_honly

      real*8  dt, aye, temstart, temend, gamma,
     &        utim, uxyz, uaye, urho, utem, fh, dtoh, z_solar, 
     &        dtemstart, dtemend, clEleFra

!  Density, energy and velocity fields fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        d(in,jn,kn),     e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &        metal(in,jn,kn), Vheat(in,jn,kn), Mheat(in,jn,kn)

!  Radiative transfer fields

      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn), photogamma(in,jn,kn)

!  Cooling tables (coolings rates as a function of temperature)

      real*8  hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec),
     &        gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), h2ltea(nratec),
     &        gasgra(nratec), ciecoa(nratec),
     &        ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec),
     &        compa, piHI, piHeI, piHeII, comp_xraya, comp_temp,
     &        gammaha 

      real*8  avgsighi, avgsighei, avgsigheii

!  Chemistry tables (rates as a function of temperature)

      real*8 k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &       k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &       k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &       k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &       k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &       k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &       k54a(nratec), k55a(nratec), k56a(nratec),
     &       k57a(nratec), k58a(nratec),
     &       k13dda(nratec, 14), h2dusta(nratec, ndratec),
     &       ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &       k24, k25, k26, k27, k28, k29, k30, k31

!  Cloudy cooling data

      integer icmbTfloor, iClHeat
      integer*8 priGridRank, priDataSize,
     &     metGridRank, metDataSize,
     &     priGridDim(priGridRank), metGridDim(metGridRank)
      real*8 priPar1(priGridDim(1)), priPar2(priGridDim(2)), 
     &     priPar3(priGridDim(3)), priPar4(priGridDim(4)),
     &     priPar5(priGridDim(5)),
     &     metPar1(metGridDim(1)), metPar2(metGridDim(2)), 
     &     metPar3(metGridDim(3)), metPar4(metGridDim(4)),
     &     metPar5(metGridDim(5)),
     &     priCooling(priDataSize), priHeating(priDataSize),
     &     priMMW(priDataSize),
     &     metCooling(metDataSize), metHeating(metDataSize)

!  Parameters

      integer itmax
      parameter (itmax = 10000)

#ifdef CONFIG_BFLOAT_4
      R_PREC tolerance
      parameter (tolerance = 1.0e-05_RKIND)
#endif

#ifdef CONFIG_BFLOAT_8
      R_PREC tolerance
      parameter (tolerance = 1.0e-10_RKIND)
#endif

      real*8 mh
      parameter (mh = mass_h)

!  Locals

      integer i, j, k, iter
      integer t, dj, dk
      real*8 ttmin, dom, energy, comp1, comp2
      real*8 coolunit, dbase1, tbase1, xbase1, chunit, uvel
      real*8 heq1, heq2, eqk221, eqk222, eqk131, eqk132,
     &       eqt1, eqt2, eqtdef, dheq, heq, dlogtem

!  row temporaries

      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in),
     &       dtit(in), ttot(in), p2d(in), tgas(in), tgasold(in),
     &       tdust(in), metallicity(in), rhoH(in), olddtit

!  Rate equation row temporaries

      real*8 HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &       HMp(in), H2Ip(in), H2IIp(in),
     &       dep(in), dedot(in),HIdot(in), dedot_prev(in),
     &       DIp(in), DIIp(in), HDIp(in), HIdot_prev(in),
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in),
     &       k13dd(in, 14), h2dust(in),
     &       ncrn(in), ncrd1(in), ncrd2(in)

!  Cooling/heating row locals

      real*8 ceHI(in), ceHeI(in), ceHeII(in),
     &       ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &       reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &       brem(in), edot(in)
      real*8 hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &       gpldl(in), gphdl(in), hdlte(in), hdlow(in), cieco(in)

!  Iteration mask

      logical itmask(in)
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set error indicator

      ierr = 0
      
!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      uvel     = uxyz / utim
c      chunit = (7.17775e-12_DKIND)/(2._DKIND*uvel*uvel*mh)   ! 4.5 eV per H2 formed
      chunit = (1.60218e-12_DKIND)/(2._DKIND*uvel*uvel*mh)   ! 1 eV per H2 formed

      dlogtem = (log(temend) - log(temstart))/real(nratec-1, DKIND)

!  Convert densities from comoving to proper

      if (iexpand .eq. 1) then

         call scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**(-3))

      endif

      call ceiling_species_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal,
     &                     is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal)


!  Loop over zones, and do an entire i-column in one go
      dk = ke - ks + 1
      dj = je - js + 1

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
! ierr is declared as shared and should be modified with atomic operation
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k, iter,
!$omp&   ttmin, energy, comp1, comp2,
!$omp&   heq1, heq2, eqk221, eqk222, eqk131, eqk132,
!$omp&   eqt1, eqt2, eqtdef, dheq, heq,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   dtit, ttot, p2d, tgas, tgasold,
!$omp&   tdust, metallicity, rhoH, olddtit,
!$omp&   HIp, HIIp, HeIp, HeIIp, HeIIIp,
!$omp&   HMp, H2Ip, H2IIp,
!$omp&   dep, dedot,HIdot, dedot_prev,
!$omp&   DIp, DIIp, HDIp, HIdot_prev,
!$omp&   k24shield, k25shield, k26shield,
!$omp&   k1 , k2 , k3 , k4 , k5,
!$omp&   k6 , k7 , k8 , k9 , k10,
!$omp&   k11, k12, k13, k14, k15,
!$omp&   k16, k17, k18, k19, k22,
!$omp&   k50, k51, k52, k53, k54,
!$omp&   k55, k56, k13dd, h2dust,
!$omp&   ncrn, ncrd1, ncrd2,
!$omp&   ceHI, ceHeI, ceHeII,
!$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
!$omp&   reHII, reHeII1, reHeII2, reHeIII,
!$omp&   brem, edot,
!$omp&   hyd01k, h2k01, vibh, roth, rotl,
!$omp&   gpldl, gphdl, hdlte, hdlow, cieco,
!$omp&   itmask )
#endif // #ifdef _OPENMP
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

!       tolerance = 1.0e-06_DKIND * dt

!     Initialize iteration mask to true for all cells.

         do i = is+1, ie+1
            itmask(i) = .true.
         enddo

!      If we are using coupled radiation with intermediate stepping,
!      set iteration mask to include only cells with radiation in the
!      intermediate coupled chemistry / energy step 
         if (iradtrans .eq. 1) then
            if (iradcoupled .eq. 1 .and. iradstep .eq. 1) then
               do i = is+1, ie+1
                  if (kphHI(i,j,k) .gt. 0) then
                      itmask(i) = .true.
                  else
                      itmask(i) = .false.
                  endif
               enddo
            endif

!      Normal rate solver, but don't double count cells with radiation
            if (iradcoupled .eq. 1 .and. iradstep .eq. 0) then
               do i = is+1, ie + 1
                  if (kphHI(i,j,k) .gt. 0) then
                     itmask(i) = .false.
                  else
                     itmask(i) = .true.
                  endif
               enddo
            endif
         endif ! end rad trans check (divergent from original code)
          


!        Set time elapsed to zero for each cell in 1D section

         do i = is+1, ie+1
            ttot(i) = 0._DKIND
         enddo

!        ------------------ Loop over subcycles ----------------

         do iter = 1, itmax

            do i = is+1, ie+1
               if (itmask(i)) then
                  dtit(i) = huge8
               endif
            enddo

!           Compute the cooling rate, tgas, tdust, and metallicity for this row

            call cool1d_multi_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool, idust, idim,
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma, fh,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                h2ltea, gasgra,
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH, 
     &                iradshield, avgsighi, avgsighei, avgsigheii,
     &                k24, k26,
     &                iradtrans, photogamma,
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat, clEleFra,
     &                priGridRank, priGridDim,
     &                priPar1, priPar2, priPar3, priPar4, priPar5,
     &                priDataSize, priCooling, priHeating, priMMW,
     &                metGridRank, metGridDim,
     &                metPar1, metPar2, metPar3, metPar4, metPar5,
     &                metDataSize, metCooling, metHeating, clnew,
     &                iVheat, iMheat, Vheat, Mheat,
     &                itmask)

            if (ispecies .gt. 0) then

!        Look-up rates as a function of temperature for 1D set of zones
!         (maybe should add itmask to this call)

            call lookup_cool_rates1d_g(temstart, temend, nratec, j, k,
     &               is, ie, ithreebody,
     &               in, jn, kn, ispecies, idust, iradshield,
     &               tgas, HI, HII, HeI, HeII, HeIII, tdust, metallicity, 
     &               k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &               k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &               k17a, k18a, k19a, k22a,
     &               k50a, k51a, k52a, k53a, k54a, k55a, k56a, 
     &               k57a, k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &               ncrna, ncrd1a, ncrd2a, 
     &               avgsighi, avgsighei, avgsigheii, piHI, piHeI,
     &               k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &               k11, k12, k13, k14, k15, k16, k17, k18,
     &               k19, k22, k24, k25, k26, 
     &               k50, k51, k52, k53, k54, k55, k56, k57,
     &               k58, k13dd, k24shield, k25shield, k26shield,
     &               h2dust, ncrn, ncrd1, ncrd2, 
     &               t1, t2, tdef, logtem, indixe, 
     &               dom, coolunit, tbase1, itmask)

!           Compute dedot and HIdot, the rates of change of de and HI
!             (should add itmask to this call)

            call rate_timestep_g(dedot, HIdot, ispecies, idust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58, 
     &                     h2dust, ncrn, ncrd1, ncrd2, rhoH, 
     &                     k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                     itmask, edot, chunit, dom)

!           Find timestep that keeps relative chemical changes below 10%

            do i = is+1, ie+1
               if (itmask(i)) then
!              Bound from below to prevent numerical errors
               
	       if (abs(dedot(i)) .lt. tiny8) 
     &             dedot(i) = min(tiny,de(i,j,k))
	       if (abs(HIdot(i)) .lt. tiny8)
     &             HIdot(i) = min(tiny,HI(i,j,k))

!              If the net rate is almost perfectly balanced then set
!                  it to zero (since it is zero to available precision)

               if (min(abs(k1(i)* de(i,j,k)*HI(i,j,k)),
     &                 abs(k2(i)*HII(i,j,k)*de(i,j,k)))/
     &             max(abs(dedot(i)),abs(HIdot(i))) .gt.
     &              1.0e6_DKIND) then
                  dedot(i) = tiny8
                  HIdot(i) = tiny8
               endif

!              If the iteration count is high then take the smaller of
!                the calculated dedot and last time step's actual dedot.
!                This is intended to get around the problem of a low
!                electron or HI fraction which is in equilibrium with high
!                individual terms (which all nearly cancel).

               if (iter .gt. 50) then
                  dedot(i) = min(abs(dedot(i)), abs(dedot_prev(i)))
                  HIdot(i) = min(abs(HIdot(i)), abs(HIdot_prev(i)))
               endif

!              compute minimum rate timestep

               olddtit = dtit(i)
               dtit(i) = min(abs(0.1_DKIND*de(i,j,k)/dedot(i)), 
     &              abs(0.1_DKIND*HI(i,j,k)/HIdot(i)),
     &              dt-ttot(i), 0.5_DKIND*dt)

               if (d(i,j,k)*dom .gt. 1e8_DKIND .and. 
     &              edot(i) .gt. 0._DKIND) then
!              Equilibrium value for H is:
!              H = (-1._DKIND / (4*k22)) * (k13 - sqrt(8 k13 k22 rho + k13^2))
!              We now want this to change by 10% or less, but we're only
!              differentiating by dT.  We have de/dt.  We need dT/de.
!              T = (g-1)*p2d*utem/N; tgas == (g-1)(p2d*utem/N)
!              dH_eq / dt = (dH_eq/dT) * (dT/de) * (de/dt)
!              dH_eq / dT (see above; we can calculate the derivative here)
!              dT / de = utem * (gamma - 1._DKIND) / N == tgas / p2d
!              de / dt = edot
!              Now we use our estimate of dT/de to get the estimated
!              difference in the equilibrium
                eqt2 = min(log(tgas(i)) + 0.1_DKIND*dlogtem, t2(i))
                eqtdef = (eqt2 - t1(i))/(t2(i) - t1(i))
                eqk222 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk132 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq2 = (-1._DKIND / (4._DKIND*eqk222)) * (eqk132-
     &               sqrt(8._DKIND*eqk132*eqk222*
     &               fh*d(i,j,k)+eqk132**2._DKIND))

                eqt1 = max(log(tgas(i)) - 0.1_DKIND*dlogtem, t1(i))
                eqtdef = (eqt1 - t1(i))/(t2(i) - t1(i))
                eqk221 = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*eqtdef
                eqk131 = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*eqtdef
                heq1 = (-1._DKIND / (4._DKIND*eqk221)) * (eqk131-
     &               sqrt(8._DKIND*eqk131*eqk221*
     &               fh*d(i,j,k)+eqk131**2._DKIND))

                dheq = (abs(heq2-heq1)/(exp(eqt2) - exp(eqt1)))
     &               * (tgas(i)/p2d(i)) * edot(i)
                heq = (-1._DKIND / (4._DKIND*k22(i))) * (k13(i)-
     &               sqrt(8._DKIND*k13(i)*k22(i)*
     &               fh*d(i,j,k)+k13(i)**2._DKIND))
                !write(0,*) heq2, heq1, eqt2, eqt1, tgas(i), p2d(i), 
!     &                     edot(i)
                if (d(i,j,k)*dom.gt.1e18_DKIND.and.i.eq.4) then
#ifdef _OPENMP
!$omp critical
                  write(0,*) HI(i,j,k)/heq, edot(i), tgas(i)
!$omp end critical
#endif
                endif
                dtit(i) = min(dtit(i), 0.1_DKIND*heq/dheq)
              endif
              if (iter.gt.10_DKIND) then
                 dtit(i) = min(olddtit*1.5_DKIND, dtit(i))
              endif

#define DONT_WRITE_COOLING_DEBUG
#ifdef WRITE_COOLING_DEBUG
!              Output some debugging information if required

               if (dtit(i)/dt .lt. 1.0e-2_DKIND .and.
     &             iter .gt. 800 .and.
     &             abs((dt-ttot(i))/dt) .gt. 1.0e-3_DKIND) then
#ifdef _OPENMP
!$omp critical
                  write(4,1000) iter,i,j,k,dtit(i),
     &              ttot(i),dt,de(i,j,k),dedot(i),HI(i,j,k),HIdot(i),
     &              tgas(i), dedot_prev(i), HIdot_prev(i)
                  write(4,1100) HI(i,j,k),HII(i,j,k),
     &              HeI(i,j,k),HeII(i,j,k),HeIII(i,j,k),
     &              HM(i,j,k),H2I(i,j,k),H2II(i,j,k),de(i,j,k)
                  write(4,1100)
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  ,
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k),
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k),
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k),
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._DKIND,
     &               - 2._DKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k),
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) ,
     &               + 2._DKIND*k13(i)*HI(i,j,k)*H2I(i,j,k)/2._DKIND,
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._DKIND,
     &               + 2._DKIND*k12(i)*de(i,j,k)*H2I(i,j,k)/2._DKIND,
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k),
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k),
     &               + 2._DKIND*k16(i)*HM(i,j,k)    *HII(i,j,k),
     &               + 2._DKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._DKIND,
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._DKIND,
     &               -      k57(i)*HI(i,j,k)    *HI(i,j,k),
     &               -      k58(i)*HI(i,j,k)    *HeI(i,j,k)/4._DKIND
!$omp end critical
#endif
               endif
 1000          format(i5,3(i3,1x),1p,11(e11.3))
 1100          format(1p,22(e11.3))
#endif /* WRITE_COOLING_DEBUG */
            else               ! itmask
               dtit(i) = dt;
            endif
            enddo               ! end loop over i

            endif ! end if (ispecies .gt. 0)

!           Compute maximum timestep for cooling/heating

            do i = is+1, ie+1
               if (itmask(i)) then
!              Set energy per unit volume of this cell based in the pressure
!              (the gamma used here is the right one even for H2 since p2d 
!               is calculated with this gamma).

               energy = max(p2d(i)/(gamma-1._DKIND), tiny8)

!              If the temperature is at the bottom of the temperature look-up 
!              table and edot < 0, then shut off the cooling.

               if (tgas(i) .le. 1.01_DKIND*temstart .and.
     &              edot(i) .lt. 0._DKIND) 
     &              edot(i) = tiny8
	       if (abs(edot(i)) .lt. tiny8) edot(i) = tiny8

!              Compute timestep for 10% change

                  dtit(i) = min(real(abs(0.1_DKIND*
     &              energy/edot(i)), DKIND), 
     &              dt-ttot(i), dtit(i))

               if (dtit(i) .ne. dtit(i)) then
#ifdef _OPENMP
!$omp critical
                 write(6,*) 'HUGE dtit :: ', energy, edot(i), dtit(i),
     &                      dt, ttot(i), abs(0.1_DKIND*energy/edot(i)), 
     &                      real(abs(0.1_DKIND*energy/edot(i)), DKIND)
!$omp end critical
#endif
               endif

#ifdef WRITE_COOLING_DEBUG
!              If the timestep is too small, then output some debugging info

               if (((dtit(i)/dt .lt. 1.0e-2_DKIND .and. iter .gt. 800) 
     &               .or. iter .gt. itmax-100) .and.
     &              abs((dt-ttot(i))/dt) .gt. 1.0e-3_DKIND) then
#ifdef _OPENMP
!$omp critical
                 write(3,2000) i,j,k,iter,e(i,j,k),edot(i),tgas(i),
     &             energy,de(i,j,k),ttot(i),d(i,j,k),e(i,j,k),dtit(i)
!$omp end critical
#endif
               endif
 2000          format(4(i4,1x),1p,10(e14.3))
#endif /* WRITE_COOLING_DEBUG */
            endif   ! itmask
            enddo   ! end loop over i

!           Update total and gas energy

            if (icool .eq. 1) then
            do i = is+1, ie+1
               if (itmask(i)) then
               e(i,j,k)  = e(i,j,k) +
     &                 real(edot(i)/d(i,j,k)*dtit(i), RKIND)
#ifdef WRITE_COOLING_DEBUG
               if (e(i,j,k) .ne. e(i,j,k)) then
#ifdef _OPENMP
!$omp critical
                    write(3,*) edot(i),d(i,j,k),dtit(i)
!$omp end critical
#endif
               endif
#endif /* WRITE_COOLING_DEBUG */

            endif               ! itmask
            enddo
            endif

            if (ispecies .gt. 0) then

!           Solve rate equations with one linearly implicit Gauss-Seidel 
!           sweep of a backward Euler method ---

            call step_rate_g(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies, idust,
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58, 
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                     itmask)

            endif

!           Add the timestep to the elapsed time for each cell and find
!            minimum elapsed time step in this row

            ttmin = huge8
            do i = is+1, ie+1
               ttot(i) = min(ttot(i) + dtit(i), dt)
               if (abs(dt-ttot(i)) .lt.
     &              0.001_DKIND*dt) itmask(i) = .false.
               if (ttot(i).lt.ttmin) ttmin = ttot(i)
            enddo

!           If all cells are done (on this slice), then exit

            if (abs(dt-ttmin) .lt. tolerance*dt) go to 9999

!           Next subcycle iteration

         enddo

 9999    continue

!       Abort if iteration count exceeds maximum

         if (iter .gt. itmax) then
#ifdef _OPENMP
!$omp critical
	    write(0,*) 'inside if statement solve rate cool:',is,ie
            write(6,*) 'MULTI_COOL iter > ',itmax,' at j,k =',j,k
            write(0,*) 'FATAL error (2) in MULTI_COOL'
            write(0,'(" dt = ",1pe10.3," ttmin = ",1pe10.3)') dt, ttmin
            write(0,'((16(1pe8.1)))') (dtit(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (ttot(i),i=is+1,ie+1)
            write(0,'((16(1pe8.1)))') (edot(i),i=is+1,ie+1)
            write(0,'((16(l3)))') (itmask(i),i=is+1,ie+1)
!$omp end critical
#endif
c            WARNING_MESSAGE
         endif

         if (iter .gt. itmax/2) then
#ifdef _OPENMP
!$omp critical
            write(6,*) 'MULTI_COOL iter,j,k =',iter,j,k
!$omp end critical
#endif
         end if
!     
!     Next j,k
!     
      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

!     Convert densities back to comoving from proper

      if (iexpand .eq. 1) then

         call scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                  HM, H2I, H2II, DI, DII, HDI, metal,
     &                  is, ie, js, je, ks, ke,
     &                  in, jn, kn, ispecies, imetal, aye**3)

      endif

      if (ispecies .gt. 0) then

!     Correct the species to ensure consistency (i.e. type conservation)

      call make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                     HM, H2I, H2II, DI, DII, HDI, metal, 
     &                     d, is, ie, js, je, ks, ke,
     &                     in, jn, kn, ispecies, imetal, fh, dtoh)

      endif

      return
      end

c -----------------------------------------------------------
!   This routine scales the density fields from comoving to
!     proper densities (and back again).

      subroutine scale_fields_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, factor)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        d(in,jn,kn),metal(in,jn,kn)
      real*8 factor

!     locals

      integer i, j, k

!     Multiply all fields by factor (1/a^3 or a^3)

      do k = ks+1, ke+1
         do j = js+1, je+1
            do i = is+1, ie+1
               d(i,j,k)     = d(i,j,k)*factor
            enddo
         enddo
      enddo

      if (ispecies .gt. 0) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  de(i,j,k)    = de(i,j,k)*factor
                  HI(i,j,k)    = HI(i,j,k)*factor
                  HII(i,j,k)   = HII(i,j,k)*factor
                  HeI(i,j,k)   = HeI(i,j,k)*factor
                  HeII(i,j,k)  = HeII(i,j,k)*factor
                  HeIII(i,j,k) = HeIII(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = HM(i,j,k)*factor
                  H2I(i,j,k)  = H2I(i,j,k)*factor
                  H2II(i,j,k) = H2II(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = DI(i,j,k)*factor
                  DII(i,j,k) = DII(i,j,k)*factor
                  HDI(i,j,k) = HDI(i,j,k)*factor
               enddo
            enddo
         enddo
      endif
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*factor
               enddo
            enddo
         enddo
      endif

      return
      end

c -----------------------------------------------------------
!   This routine ensures that the species aren't below tiny.

      subroutine ceiling_species_g(d, de, HI, HII, HeI, HeII, HeIII,
     &                           HM, H2I, H2II, DI, DII, HDI, metal,
     &                           is, ie, js, je, ks, ke,
     &                           in, jn, kn, ispecies, imetal)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC  d(in,jn,kn),
     &        de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)

!     locals

      integer i, j, k

      if (ispecies .gt. 0) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  de(i,j,k)    = max(de(i,j,k), tiny)
                  HI(i,j,k)    = max(HI(i,j,k), tiny)
                  HII(i,j,k)   = max(HII(i,j,k), tiny)
                  HeI(i,j,k)   = max(HeI(i,j,k), tiny)
                  HeII(i,j,k)  = max(HeII(i,j,k), tiny)
                  HeIII(i,j,k) = max(HeIII(i,j,k), 1e-5_RKIND*tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  HM(i,j,k)   = max(HM(i,j,k), tiny)
                  H2I(i,j,k)  = max(H2I(i,j,k), tiny)
                  H2II(i,j,k) = max(H2II(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (ispecies .gt. 2) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  DI(i,j,k)  = max(DI(i,j,k), tiny)
                  DII(i,j,k) = max(DII(i,j,k), tiny)
                  HDI(i,j,k) = max(HDI(i,j,k), tiny)
               enddo
            enddo
         enddo
      endif
      if (imetal .eq. 1) then
         do k = ks+1, ke+1
            do j = js+1, je+1
               do i = is+1, ie+1
                  metal(i,j,k) = min(max(metal(i,j,k), tiny),
     &                 0.9_RKIND*d(i,j,k))
               enddo
            enddo
         enddo
      endif

      return
      end



! -----------------------------------------------------------
! This routine uses the temperature to look up the chemical
!   rates which are tabulated in a log table as a function
!   of temperature.

      subroutine lookup_cool_rates1d_g(temstart, temend, nratec, j, k,
     &                is, ie, ithreebody,
     &                in, jn, kn, ispecies, idust, iradshield,
     &                tgas1d, HI, HII, HeI, HeII, HeIII, tdust, metallicity, 
     &                k1a, k2a, k3a, k4a, k5a, k6a, k7a, k8a, k9a, k10a,
     &                k11a, k12a, k13a, k13dda, k14a, k15a, k16a,
     &                k17a, k18a, k19a, k22a,
     &                k50a, k51a, k52a, k53a, k54a, k55a, k56a, k57a,
     &                k58a, ndratec, dtemstart, dtemend, h2dusta, 
     &                ncrna, ncrd1a, ncrd2a,
     &                avgsighi, avgsighei, avgsigheii, piHI, piHeI,
     &                k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     &                k11, k12, k13, k14, k15, k16, k17, k18,
     &                k19, k22, k24, k25, k26,
     &                k50, k51, k52, k53, k54, k55, k56, k57,
     &                k58, k13dd, k24shield, k25shield, k26shield,
     &                h2dust, ncrn, ncrd1, ncrd2,
     &                t1, t2, tdef, logtem, indixe, 
     &                dom, coolunit, tbase1, itmask)
! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer is, ie, nratec,
     &        in, jn, kn, ispecies, ithreebody, j, k,
     &        idust, ndratec, iradshield
      real*8 temstart, temend, tgas1d(in), dom,
     &       dtemstart, dtemend
      real*8 coolunit, tbase1
      logical itmask(in)

!     Chemistry rates as a function of temperature

      real*8 k1a (nratec), k2a (nratec), k3a (nratec), k4a (nratec), 
     &       k5a (nratec), k6a (nratec), k7a (nratec), k8a (nratec), 
     &       k9a (nratec), k10a(nratec), k11a(nratec), k12a(nratec), 
     &       k13a(nratec), k14a(nratec), k15a(nratec), k16a(nratec), 
     &       k17a(nratec), k18a(nratec), k19a(nratec), k22a(nratec),
     &       k50a(nratec), k51a(nratec), k52a(nratec), k53a(nratec),
     &       k54a(nratec), k55a(nratec), k56a(nratec), k57a(nratec),
     &       k58a(nratec), k13dda(nratec, 14), h2dusta(nratec, ndratec), 
     &       ncrna(nratec), ncrd1a(nratec), ncrd2a(nratec),
     &       k24, k25, k26,
     &       piHI, piHeI,
     &       avgsighi, avgsighei, avgsigheii

!     Density fields

      R_PREC  HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)

!     Returned rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in),
     &       k13dd(in, 14), h2dust(in), 
     &       ncrn(in), ncrd1(in), ncrd2(in),
     &       k24shield(in), k25shield(in), k26shield(in)

!     1D temporaries (passed in)

      integer*8 indixe(in)
      real*8 t1(in), t2(in), logtem(in), tdef(in),
     &       tdust(in), metallicity(in)

!     1D temporaries (not passed in)

      integer*8 d_indixe(in)
      real*8 d_t1(in), d_t2(in), d_logtem(in), d_tdef(in),
     &       dusti1(in), dusti2(in)

!     Parameters

      real*8 everg, e24, e26
      parameter(everg = ev2erg, e24 = 13.6_DKIND,
     &     e26 = 24.6_DKIND)

!     locals

      integer i, n1
      real*8 factor, x, logtem0, logtem9, dlogtem, nh,
     &       d_logtem0, d_logtem9, d_dlogtem
      real*8 k13_CID, k13_DT

      real*8 nSSh, nratio

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem = (log(temend) - log(temstart))/real(nratec-1, DKIND)

      do i = is+1, ie+1
         if (itmask(i)) then
!        Compute temp-centered temperature (and log)

!        logtem(i) = log(0.5_DKIND*(tgas(i)+tgasold(i)))
         logtem(i) = log(tgas1d(i))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Find index into tble and precompute interpolation values

         indixe(i) = min(nratec-1,
     &        max(1,int((logtem(i)-logtem0)/dlogtem, DIKIND)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!        Do linear table lookup (in log temperature)

         k1(i) = k1a(indixe(i)) +
     &           (k1a(indixe(i)+1) -k1a(indixe(i)))*tdef(i)
         k2(i) = k2a(indixe(i)) +
     &           (k2a(indixe(i)+1) -k2a(indixe(i)))*tdef(i)
         k3(i) = k3a(indixe(i)) +
     &           (k3a(indixe(i)+1) -k3a(indixe(i)))*tdef(i)
         k4(i) = k4a(indixe(i)) +
     &           (k4a(indixe(i)+1) -k4a(indixe(i)))*tdef(i)
         k5(i) = k5a(indixe(i)) +
     &           (k5a(indixe(i)+1) -k5a(indixe(i)))*tdef(i)
         k6(i) = k6a(indixe(i)) +
     &           (k6a(indixe(i)+1) -k6a(indixe(i)))*tdef(i)
         k57(i) = k57a(indixe(i)) +
     &            (k57a(indixe(i)+1) -k57a(indixe(i)))*tdef(i)
         k58(i) = k58a(indixe(i)) +
     &            (k58a(indixe(i)+1) -k58a(indixe(i)))*tdef(i)
      endif
      enddo

!     Look-up for 9-species model

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k7(i) = k7a(indixe(i)) +
     &            (k7a(indixe(i)+1) -k7a(indixe(i)))*tdef(i)
            k8(i) = k8a(indixe(i)) +
     &            (k8a(indixe(i)+1) -k8a(indixe(i)))*tdef(i)
            k9(i) = k9a(indixe(i)) +
     &            (k9a(indixe(i)+1) -k9a(indixe(i)))*tdef(i)
            k10(i) = k10a(indixe(i)) +
     &            (k10a(indixe(i)+1) -k10a(indixe(i)))*tdef(i)
            k11(i) = k11a(indixe(i)) +
     &            (k11a(indixe(i)+1) -k11a(indixe(i)))*tdef(i)
            k12(i) = k12a(indixe(i)) +
     &            (k12a(indixe(i)+1) -k12a(indixe(i)))*tdef(i)
            k13(i) = k13a(indixe(i)) +
     &            (k13a(indixe(i)+1) -k13a(indixe(i)))*tdef(i)
            k14(i) = k14a(indixe(i)) +
     &            (k14a(indixe(i)+1) -k14a(indixe(i)))*tdef(i)
            k15(i) = k15a(indixe(i)) +
     &            (k15a(indixe(i)+1) -k15a(indixe(i)))*tdef(i)
            k16(i) = k16a(indixe(i)) +
     &            (k16a(indixe(i)+1) -k16a(indixe(i)))*tdef(i)
            k17(i) = k17a(indixe(i)) +
     &            (k17a(indixe(i)+1) -k17a(indixe(i)))*tdef(i)
            k18(i) = k18a(indixe(i)) +
     &            (k18a(indixe(i)+1) -k18a(indixe(i)))*tdef(i)
            k19(i) = k19a(indixe(i)) +
     &            (k19a(indixe(i)+1) -k19a(indixe(i)))*tdef(i)
            k22(i) = k22a(indixe(i)) +
     &            (k22a(indixe(i)+1) -k22a(indixe(i)))*tdef(i)

!     H2 formation heating terms.

            ncrn(i) = ncrna(indixe(i)) +
     &           (ncrna(indixe(i)+1) -ncrna(indixe(i)))*tdef(i)
            ncrd1(i) = ncrd1a(indixe(i)) +
     &           (ncrd1a(indixe(i)+1) -ncrd1a(indixe(i)))*tdef(i)
            ncrd2(i) = ncrd2a(indixe(i)) +
     &           (ncrd2a(indixe(i)+1) -ncrd2a(indixe(i)))*tdef(i)

         endif
         enddo

         do n1 = 1, 14
            do i = is+1, ie+1
               if (itmask(i)) then
               k13dd(i,n1) = k13dda(indixe(i),n1) +
     &             (k13dda(indixe(i)+1,n1) - 
     &               k13dda(indixe(i)  ,n1) )*tdef(i)
            endif
            enddo
         enddo         

      endif

!     Look-up for 12-species model

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
            k50(i) = k50a(indixe(i)) +
     &            (k50a(indixe(i)+1) -k50a(indixe(i)))*tdef(i)
            k51(i) = k51a(indixe(i)) +
     &            (k51a(indixe(i)+1) -k51a(indixe(i)))*tdef(i)
            k52(i) = k52a(indixe(i)) +
     &            (k52a(indixe(i)+1) -k52a(indixe(i)))*tdef(i)
            k53(i) = k53a(indixe(i)) +
     &            (k53a(indixe(i)+1) -k53a(indixe(i)))*tdef(i)
            k54(i) = k54a(indixe(i)) +
     &            (k54a(indixe(i)+1) -k54a(indixe(i)))*tdef(i)
            k55(i) = k55a(indixe(i)) +
     &            (k55a(indixe(i)+1) -k55a(indixe(i)))*tdef(i)
            k56(i) = k56a(indixe(i)) +
     &            (k56a(indixe(i)+1) -k56a(indixe(i)))*tdef(i)
         endif
         enddo
      endif

!     Look-up for H2 formation on dust

      if (idust .gt. 0) then

         d_logtem0 = log(dtemstart)
         d_logtem9 = log(dtemend)
         d_dlogtem = (log(dtemend) - log(dtemstart))/
     &        real(ndratec-1, DKIND)

         do i = is+1, ie+1
            if (itmask(i)) then

!              Assume dust melting at T > 1500 K

               if (tdust(i) .gt. dtemend) then
                  h2dust(i) = tiny8
               else

!                 Get log dust temperature

                  d_logtem(i) = log(tdust(i))
                  d_logtem(i) = max(d_logtem(i), d_logtem0)
                  d_logtem(i) = min(d_logtem(i), d_logtem9)

!                 Find index into table and precompute interpolation values

                  d_indixe(i) = min(ndratec-1,
     &                 max(1,
     &                 int((d_logtem(i)-d_logtem0)/d_dlogtem,
     &                 DIKIND)+1))
                  d_t1(i) = (d_logtem0 + (d_indixe(i) - 1)*d_dlogtem)
                  d_t2(i) = (d_logtem0 + (d_indixe(i)    )*d_dlogtem)
                  d_tdef(i) = (d_logtem(i) - d_t1(i)) / 
     &                 (d_t2(i) - d_t1(i))

!                 Get rate from 2D interpolation

                  dusti1(i) = h2dusta(indixe(i), d_indixe(i)) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)) - 
     &                 h2dusta(indixe(i),   d_indixe(i)))*tdef(i)
                  dusti2(i) = h2dusta(indixe(i), d_indixe(i)+1) +
     &                 (h2dusta(indixe(i)+1, d_indixe(i)+1) - 
     &                 h2dusta(indixe(i),   d_indixe(i)+1))*tdef(i)
                  h2dust(i) = dusti1(i) + 
     &                 (dusti2(i) - dusti1(i))*d_tdef(i)

!                 Multiply by metallicity

                  h2dust(i) = h2dust(i) * metallicity(i)

               endif
            endif
         enddo
      endif

!        Include approximate self-shielding factors if requested

      do i = is+1, ie+1
         if (itmask(i)) then
            k24shield(i) = k24
            k25shield(i) = k25
            k26shield(i) = k26
         endif
      enddo

      if (iradshield == 1) then
!
!     approximate self shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     to shield HI, while leaving HeI and HeII optically thin
!
        do i = is+1, ie+1
          if (itmask(i)) then
            nSSh = 6.73e-3_DKIND *
     &        (avgsighi /2.49e-18_DKIND)**(-2._DKIND/3._DKIND) *
     &        (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND) *
     &        (k24/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)
            nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh
            k24shield(i)=k24shield(i)*(0.98_DKIND*
     &            (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &             0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

             k25shield(i) = k25
             k26shield(i) = k26
          endif
        enddo

      else if (iradshield == 2) then
!
!     Better self-shielding in HI using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     approximate self shielding in HeI and HeII
!
        do i = is+1, ie+1
          if (itmask(i)) then
            nSSh = 6.73e-3_DKIND *
     &        (avgsighi/2.49e-18_DKIND)**(-2._DKIND/3._DKIND) *
     &        (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND) *
     &        (k24/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)
            nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh

            k24shield(i)=k24shield(i)*(0.98_DKIND*
     &            (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &             0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

!
!     Apply same equations to HeI (assumes HeI closely follows HI)
!
            nSSh = 6.73e-3_DKIND *
     &        (avgsighei/2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &        (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &        (k26/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)
            nratio = 0.25_DKIND*
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))*dom/nSSh

            k26shield(i)=k26shield(i)*(0.98_DKIND*
     &            (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &             0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

            k25shield(i) = k25
          endif
        enddo

      else if (iradshield == 3) then
!
!     shielding using Eq. 13 and 14 from
!     Rahmati et. al. 2013 (MNRAS, 430, 2427-2445)
!     in HI and HeI, but ignoring HeII heating entirely
!
        do i = is+1, ie+1
          if (itmask(i)) then
            nSSh = 6.73e-3_DKIND *
     &        (avgsighi/2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &        (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND) *
     &        (k24/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)
            nratio = (HI(i,j,k) + HII(i,j,k))*dom/nSSh

            k24shield(i)=k24shield(i)*(0.98_DKIND*
     &            (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &             0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))

!
!     Apply same equations to HeI (assumes HeI closely follows HI)
!
            nSSh = 6.73e-3_DKIND *
     &        (avgsighei/2.49e-18_DKIND)**(-2._DKIND/3._DKIND)*
     &        (tgas1d(i)/1.0e4_DKIND)**(0.17_DKIND)*
     &        (k26/tbase1/1.0e-12_DKIND)**(2.0_DKIND/3.0_DKIND)
            nratio = 0.25_DKIND *
     &            (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))*dom/nSSh

            k26shield(i)=k26shield(i)*(0.98_DKIND*
     &            (1.0_DKIND+nratio**(1.64_DKIND))**(-2.28_DKIND) +
     &             0.02_DKIND*(1.0_DKIND+nratio)**(-0.84_DKIND))


            k25shield(i) = 0.0_DKIND

          endif
        enddo

      endif

#ifdef SECONDARY_IONIZATION_NOT_YET_IMPLEMENTED
!        If using a high-energy radiation field, then account for
!          effects of secondary electrons (Shull * Steenberg 1985)
!          (see calc_rate.src)

      do i = is+1, ie+1
         if (itmask(i)) then
            x = max(HII(i,j,k)/(HI(i,j,k)+HII(i,j,k)), 1.0e-4_DKIND)
            factor = 0.3908_DKIND*(1._DKIND -
     &           x**0.4092_DKIND)**1.7592_DKIND
            k24shield(i) = k24shield(i) + 
     &           factor*(piHI + 0.08_DKIND*piHeI)/
     &           (e24*everg) *coolunit*tbase1
            factor = 0.0554_DKIND*
     &           (1._DKIND - x**0.4614_DKIND)**1.6660_DKIND
            k26shield(i) = k26shield(i) + 
     &           factor*(piHI/0.08_DKIND + piHeI)/
     &           (e26*everg) *coolunit*tbase1
         endif
      enddo
#endif


!           If using H2, and using the density-dependent collisional
!             H2 dissociation rate, then replace the the density-independant
!                k13 rate with the new one.
!         May/00: there appears to be a problem with the density-dependent
!             collisional rates.  Currently turned off until further notice.

#define USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
#ifdef USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE
            if (ispecies .gt. 1 .and. ithreebody .eq. 0) then
               do i = is+1, ie+1
                  if (itmask(i)) then
                  nh = min(HI(i,j,k)*dom, 1.0e9_DKIND)
                  k13(i) = tiny8
                  if (tgas1d(i) .ge. 500._DKIND .and.
     &                tgas1d(i) .lt. 1.0e6_DKIND) then
c Direct collisional dissociation
                     k13_CID = k13dd(i,1)-k13dd(i,2)/
     &                          (1._DKIND+(nh/k13dd(i,5))**k13dd(i,7))
     &                     + k13dd(i,3)-k13dd(i,4)/
     &                          (1._DKIND+(nh/k13dd(i,6))**k13dd(i,7))
                     k13_CID = max(10._DKIND**k13_CID, tiny8)
c Dissociative tunnelling
                     k13_DT  = k13dd(i,8)-k13dd(i,9)/
     &                          (1._DKIND+(nh/k13dd(i,12))**k13dd(i,14))
     &                     + k13dd(i,10)-k13dd(i,11)/
     &                          (1._DKIND+(nh/k13dd(i,13))**k13dd(i,14))
                     k13_DT  = max(10._DKIND**k13_DT, tiny8)
c
                     k13(i)  = k13_DT + k13_CID
                  endif
               endif
               enddo
            endif
#endif /*  USE_DENSITY_DEPENDENT_H2_DISSOCIATION_RATE */

      return
      end

! -------------------------------------------------------------------
!  This routine calculates the electron and HI rates of change in
!    order to determine the maximum permitted timestep

      subroutine rate_timestep_g(dedot, HIdot, ispecies, idust,
     &                     de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II,
     &                     in, jn, kn, is, ie, j, k, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58, 
     &                     h2dust, ncrn, ncrd1, ncrd2, rhoH, 
     &                     k24shield, k25shield, k26shield,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                     itmask, edot, chunit, dom)

! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     arguments

      integer ispecies, idust, is, ie, j, k, in, jn, kn,
     &        iradtrans, irt_honly
      real*8 dedot(in), HIdot(in), dom
      real*8 edot(in)
      logical itmask(in)

!     Density fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn)

!      Radiative Transfer Fields
      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)

      real*8 chunit

!     Rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in), h2dust(in), 
     &       ncrn(in), ncrd1(in), ncrd2(in), rhoH(in), 
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k24, k25, k26, k27, k28, k29, k30, k31

!     locals

      integer i
      real*8 h2heatfac(in), H2delta(in), H2dmag, atten, tau

      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then
!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k3(i)*HeI(i,j,k)*de(i,j,k)/4._DKIND
     &               + k5(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND
     &               - k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k4(i)*HeII(i,j,k)*de(i,j,k)/4._DKIND
     &               - k6(i)*HeIII(i,j,k)*de(i,j,k)/4._DKIND
     &               + k57(i)*HI(i,j,k)*HI(i,j,k)
     &               + k58(i)*HI(i,j,k)*HeI(i,j,k)/4._DKIND
     &               +      ( k24shield(i)*HI(i,j,k)
     &               + k25shield(i)*HeII(i,j,k)/4._DKIND
     &               + k26shield(i)*HeI(i,j,k)/4._DKIND)

!     Compute the HI density rate-of-change

            HIdot(i) =
     &               - k1(i)*HI(i,j,k)*de(i,j,k)
     &               + k2(i)*HII(i,j,k)*de(i,j,k)
     &               - k57(i)*HI(i,j,k)*HI(i,j,k)
     &               - k58(i)*HI(i,j,k)*HeI(i,j,k)/4._DKIND
     &               -      k24shield(i)*HI(i,j,k)

         endif                  ! itmask
         enddo
      else

!         Include molecular hydrogen rates for HIdot

         do i = is+1, ie+1
            if (itmask(i)) then
               HIdot(i) = 
     &               -      k1(i) *de(i,j,k)    *HI(i,j,k)  
     &               -      k7(i) *de(i,j,k)    *HI(i,j,k)
     &               -      k8(i) *HM(i,j,k)    *HI(i,j,k)
     &               -      k9(i) *HII(i,j,k)   *HI(i,j,k)
     &               -      k10(i)*H2II(i,j,k)  *HI(i,j,k)/2._DKIND
     &               - 2._DKIND*k22(i)*HI(i,j,k)**2 *HI(i,j,k)
     &               +      k2(i) *HII(i,j,k)   *de(i,j,k) 
     &               + 2._DKIND*k13(i)*HI(i,j,k)    *H2I(i,j,k)/2._DKIND
     &               +      k11(i)*HII(i,j,k)   *H2I(i,j,k)/2._DKIND
     &               + 2._DKIND*k12(i)*de(i,j,k)    *H2I(i,j,k)/2._DKIND
     &               +      k14(i)*HM(i,j,k)    *de(i,j,k)
     &               +      k15(i)*HM(i,j,k)    *HI(i,j,k)
     &               + 2._DKIND*k16(i)*HM(i,j,k)    *HII(i,j,k)
     &               + 2._DKIND*k18(i)*H2II(i,j,k)  *de(i,j,k)/2._DKIND
     &               +      k19(i)*H2II(i,j,k)  *HM(i,j,k)/2._DKIND
     &               -      k57(i)*HI(i,j,k)    *HI(i,j,k)
     &               -      k58(i)*HI(i,j,k)    *HeI(i,j,k)/4._DKIND
     &               -      k24shield(i)*HI(i,j,k)

!     Add H2 formation on dust grains

            if (idust .gt. 0) then
               HIdot(i) = HIdot(i) 
     &              - 2._DKIND * h2dust(i) * rhoH(i)
            endif

!     Compute the electron density rate-of-change

            dedot(i) = 
     &               + k1(i) * HI(i,j,k)   * de(i,j,k)
     &               + k3(i) * HeI(i,j,k)  * de(i,j,k)/4._DKIND
     &               + k5(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND
     &               + k8(i) * HM(i,j,k)   * HI(i,j,k)
     &               + k15(i)* HM(i,j,k)   * HI(i,j,k)
     &               + k17(i)* HM(i,j,k)   * HII(i,j,k)
     &               + k14(i)* HM(i,j,k)   * de(i,j,k)
     &               - k2(i) * HII(i,j,k)  * de(i,j,k)
     &               - k4(i) * HeII(i,j,k) * de(i,j,k)/4._DKIND
     &               - k6(i) * HeIII(i,j,k)* de(i,j,k)/4._DKIND
     &               - k7(i) * HI(i,j,k)   * de(i,j,k)
     &               - k18(i)* H2II(i,j,k) * de(i,j,k)/2._DKIND
     &               + k57(i)* HI(i,j,k)   * HI(i,j,k)
     &               + k58(i)* HI(i,j,k)   * HeI(i,j,k)/4._DKIND
     &               + (k24shield(i)*HI(i,j,k)
     &               +  k25shield(i)*HeII(i,j,k)/4._DKIND
     &               +  k26shield(i)*HeI(i,j,k)/4._DKIND)

!     H2 formation heating

!     Equation 23 from Omukai (2000)
            h2heatfac(i) = (1._DKIND + (ncrn(i) / (dom *
     &           (HI(i,j,k) * ncrd1(i) +
     &           H2I(i,j,k) * 0.5_DKIND * ncrd2(i)))))**(-1._DKIND)

            H2delta(i) = 
     &          HI(i,j,k) *
     &           (  4.48_DKIND * k22(i) * HI(i,j,k)**2._DKIND
     &            - 4.48_DKIND * k13(i) * H2I(i,j,k)/2._DKIND)

            ! We only want to apply this if the formation dominates, but we
            ! need to apply it outside the delta calculation.
            if(H2delta(i).gt.0._DKIND) then
              H2delta(i) = H2delta(i) * h2heatfac(i)
            endif

            if (idust .gt. 0) then
               H2delta(i) = H2delta(i) + 
     &              h2dust(i) * HI(i,j,k) * rhoH(i) * 
     &              (0.2_DKIND + 4.2_DKIND * h2heatfac(i))
            endif

!            H2dmag = abs(H2delta)/(
!     &          HI(i,j,k)*( k22(i) * HI(i,j,k)**2._DKIND
!     &                    + k13(i) * H2I(i,j,k)/2._DKIND))
!            tau = (H2dmag/1e-5_DKIND)**-1.0_DKIND
!            tau = max(tau, 1.e-5_DKIND)
!            atten = min((1.-exp(-tau))/tau,1._DKIND)
            atten = 1._DKIND
            edot(i) = edot(i) + chunit * H2delta(i) * atten
!     &       + H2I(i,j,k)*( k21(i) * HI(i,j,k)**2.0_DKIND
!     &                    - k23(i) * H2I(i,j,k))
!H * (k22 * H^2 - k13 * H_2) + H_2 * (k21 * H^2 - k23 * H_2) */
         endif                  ! itmask
         enddo
      endif

!     Add photo-ionization rates if needed

      if (iradtrans .eq. 1) then
         if (irt_honly .eq. 0) then
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
     &                 + kphHeI(i,j,k) * HeI(i,j,k) / 4._DKIND
     &                 + kphHeII(i,j,k) * HeII(i,j,k) / 4._DKIND
               endif
            enddo
         else
            do i = is+1, ie+1
               if (itmask(i)) then
                  HIdot(i) = HIdot(i) - kphHI(i,j,k)*HI(i,j,k)
                  dedot(i) = dedot(i) + kphHI(i,j,k)*HI(i,j,k)
               endif
            enddo
         endif
      endif

      


      return
      end


! -----------------------------------------------------------
!  This routine uses one linearly implicit Gauss-Seidel sweep of 
!   a backward-Euler time integrator to advance the rate equations 
!   by one (sub-)cycle (dtit).

      subroutine step_rate_g(de, HI, HII, HeI, HeII, HeIII, d,
     &                     HM, H2I, H2II, DI, DII, HDI, dtit,
     &                     in, jn, kn, is, ie, j, k, ispecies, idust, 
     &                     k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11,
     &                     k12, k13, k14, k15, k16, k17, k18, k19, k22,
     &                     k24, k25, k26, k27, k28, k29, k30, k31,
     &                     k50, k51, k52, k53, k54, k55, k56, k57, k58,
     &                     h2dust, rhoH,
     &                     k24shield, k25shield, k26shield,
     &                     HIp, HIIp, HeIp, HeIIp, HeIIIp, dep,
     &                     HMp, H2Ip, H2IIp, DIp, DIIp, HDIp,
     &                     dedot_prev, HIdot_prev,
     &                     iradtrans, irt_honly, kphHI, kphHeI, kphHeII, kdissH2I,
     &                     itmask)
c -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     arguments

      integer ispecies, idust, in, jn, kn, is, ie, j, k,
     &        iradtrans, irt_honly
      real*8 dtit(in), dedot_prev(in), HIdot_prev(in)
      logical itmask(in)

!     Density fields

      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)

!     Radiation Fields
      R_PREC  kphHI(in,jn,kn), kphHeI(in,jn,kn), kphHeII(in,jn,kn),
     &        kdissH2I(in,jn,kn)


!     Rate values

      real*8 k1 (in), k2 (in), k3 (in), k4 (in), k5 (in),
     &       k6 (in), k7 (in), k8 (in), k9 (in), k10(in),
     &       k11(in), k12(in), k13(in), k14(in), k15(in),
     &       k16(in), k17(in), k18(in), k19(in), k22(in),
     &       k50(in), k51(in), k52(in), k53(in), k54(in),
     &       k55(in), k56(in), k57(in), k58(in), 
     &       h2dust(in), rhoH(in),
     &       k24shield(in), k25shield(in), k26shield(in),
     &       k24, k25, k26, k27, k28, k29, k30, k31

!     temporaries (passed in)

      real*8 HIp(in), HIIp(in), HeIp(in), HeIIp(in), HeIIIp(in),
     &       HMp(in), H2Ip(in), H2IIp(in), dep(in),
     &       DIp(in), DIIp(in), HDIp(in)

!     locals

      integer i
      real*8 scoef, acoef

!   A) the 6-species integrator
!      
      if (ispecies .eq. 1) then

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI

            scoef  = k2(i)*HII(i,j,k)*de(i,j,k)
            acoef  = k1(i)*de(i,j,k)
     &             + k57(i)*HI(i,j,k)
     &             + k58(i)*HeI(i,j,k)/4._DKIND
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            HIp(i)  = (scoef*dtit(i) + HI(i,j,k))/
     &           (1._DKIND + acoef*dtit(i))
            if (HIp(i) .ne. HIp(i)) then
#ifdef _OPENMP
!$omp critical
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), kphHI(i,j,k)
     $              scoef, acoef, dtit(i)
!$omp end critical
#endif
c               ERROR_MESSAGE
            endif

!        2) HII
c 
            scoef  = k1(i)*HIp(i)*de(i,j,k)
     &             + k57(i)*HIp(i)*HIp(i)
     &             + k58(i)*HIp(i)*HeI(i,j,k)/4._DKIND
     &             + k24shield(i)*HIp(i)
            if (iradtrans .eq. 1) 
     &          scoef = scoef + kphHI(i,j,k)*HIp(i)
            acoef  = k2(i)*de (i,j,k)
            HIIp(i) = (scoef*dtit(i) + HII(i,j,k))/
     &           (1._DKIND +acoef*dtit(i))
!
            if (HIIp(i) .le. 0._DKIND) then  !#####
#ifdef _OPENMP
!$omp critical
               write(*,*) 'negative HIIp! :: ', i, j, k, HIIp(i), 
     $              scoef, dtit(i), HII(i,j,k), acoef,
     $              k2(i), de(i,j,k),
     $              kphHI(i,j,k), HIp(i),
     $              k24shield(i)
!$omp end critical
#endif
            endif

!        3) Electron density

            scoef = 0._DKIND
     &                 + k57(i)*HIp(i)*HIp(i)
     &                 + k58(i)*HIp(i)*HeI(i,j,k)/4._DKIND
     &                 + k24shield(i)*HI(i,j,k)
     &                 + k25shield(i)*HeII(i,j,k)/4._DKIND
     &                 + k26shield(i)*HeI(i,j,k)/4._DKIND

            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 0) )
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)
     &                + kphHeI(i,j,k)  * HeI(i,j,k)  / 4._DKIND
     &                + kphHeII(i,j,k) * HeII(i,j,k) / 4._DKIND
            if ( (iradtrans .eq. 1) .and. ( irt_honly .eq. 1) )
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)



            acoef = -(k1(i)*HI(i,j,k)      - k2(i)*HII(i,j,k)
     &              + k3(i)*HeI(i,j,k)/4._DKIND -
     &           k6(i)*HeIII(i,j,k)/4._DKIND
     &              + k5(i)*HeII(i,j,k)/4._DKIND -
     &           k4(i)*HeII(i,j,k)/4._DKIND)
            dep(i)   = (scoef*dtit(i) + de(i,j,k))
     &                     / (1._DKIND + acoef*dtit(i))

         endif                  ! itmask
         enddo

      endif                     ! (ispecies .eq. 1)

!  --- (B) Do helium chemistry in any case: (for all ispecies values) ---

      do i = is+1, ie+1
         if (itmask(i)) then

!        4) HeI

         scoef  = k4(i)*HeII(i,j,k)*de(i,j,k)
         acoef  = k3(i)*de(i,j,k)
     &                + k26shield(i)

         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeI(i,j,k)

         HeIp(i)   = ( scoef*dtit(i) + HeI(i,j,k) ) 
     &              / ( 1._DKIND + acoef*dtit(i) )

!        5) HeII

         scoef  = k3(i)*HeIp(i)*de(i,j,k)
     &          + k6(i)*HeIII(i,j,k)*de(i,j,k)
     &          + k26shield(i)*HeIp(i)
     
         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeI(i,j,k)*HeIp(i)

         acoef  = k4(i)*de(i,j,k) + k5(i)*de(i,j,k)
     &          + k25shield(i)
     
         if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       acoef = acoef + kphHeII(i,j,k)

         HeIIp(i)  = ( scoef*dtit(i) + HeII(i,j,k) )
     &              / ( 1._DKIND + acoef*dtit(i) )

!       6) HeIII

         scoef   = k5(i)*HeIIp(i)*de(i,j,k)
     &           + k25shield(i)*HeIIp(i)
         if ((iradtrans .eq. 1) .and. (irt_honly .eq. 0))
     &       scoef = scoef + kphHeII(i,j,k) * HeIIp(i)
         acoef   = k6(i)*de(i,j,k)
         HeIIIp(i)  = ( scoef*dtit(i) + HeIII(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

      endif                     ! itmask
      enddo

c --- (C) Now do extra 3-species for molecular hydrogen ---

      if (ispecies .gt. 1) then

!        First, do HI/HII with molecular hydrogen terms

         do i = is+1, ie+1
            if (itmask(i)) then

!        1) HI
!     
            scoef  =      k2(i) * HII(i,j,k) * de(i,j,k) 
     &             + 2._DKIND*k13(i)* HI(i,j,k)  * H2I(i,j,k)/2._DKIND
     &             +      k11(i)* HII(i,j,k) * H2I(i,j,k)/2._DKIND
     &             + 2._DKIND*k12(i)* de(i,j,k)  * H2I(i,j,k)/2._DKIND
     &             +      k14(i)* HM(i,j,k)  * de(i,j,k)
     &             +      k15(i)* HM(i,j,k)  * HI(i,j,k)
     &             + 2._DKIND*k16(i)* HM(i,j,k)  * HII(i,j,k)
     &             + 2._DKIND*k18(i)* H2II(i,j,k)* de(i,j,k)/2._DKIND
     &             +      k19(i)* H2II(i,j,k)* HM(i,j,k)/2._DKIND
     &             + 2._DKIND*k31   * H2I(i,j,k)/2._DKIND
            if (iradtrans .eq. 1)
     &      scoef = scoef+2._DKIND* kdissH2I(i,j,k)* H2I(i,j,k)/2._DKIND
            acoef  =      k1(i) * de(i,j,k)
     &             +      k7(i) * de(i,j,k)  
     &             +      k8(i) * HM(i,j,k)
     &             +      k9(i) * HII(i,j,k)
     &             +      k10(i)* H2II(i,j,k)/2._DKIND
     &             + 2._DKIND*k22(i)* HI(i,j,k)**2
     &             +      k57(i)* HI(i,j,k)
     &             +      k58(i)* HeI(i,j,k)/4._DKIND
     &             + k24shield(i)

            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)

            if (idust .gt. 0) then
               acoef = acoef + 2._DKIND * h2dust(i) * rhoH(i)
            endif

            HIp(i)  = ( scoef*dtit(i) + HI(i,j,k) ) / 
     &                      ( 1. + acoef*dtit(i) )
            if (HIp(i) .ne. HIp(i)) then
#ifdef _OPENMP
!$omp critical
               write(*,*) 'HUGE HIp! :: ', i, j, k, HIp(i), HI(i,j,k),
     $              HII(i,j,k), de(i,j,k), H2I(i,j,k), kdissH2I(i,j,k),
     $              kphHI(i,j,k)
!$omp end critical
#endif
            endif

!          2) HII

            scoef  =    k1(i)  * HI(i,j,k) * de(i,j,k)
     &             +    k10(i) * H2II(i,j,k)*HI(i,j,k)/2._DKIND
     &             +    k57(i) * HI(i,j,k) * HI(i,j,k)
     &             +    k58(i) * HI(i,j,k) * HeI(i,j,k)/4._DKIND
     &             + k24shield(i)*HI(i,j,k)
            if (iradtrans .eq. 1) 
     &          scoef = scoef + kphHI(i,j,k) * HI(i,j,k)

            acoef  =    k2(i)  * de(i,j,k)
     &             +    k9(i)  * HI(i,j,k)
     &             +    k11(i) * H2I(i,j,k)/2._DKIND
     &             +    k16(i) * HM(i,j,k)
     &             +    k17(i) * HM(i,j,k)
            HIIp(i)   = ( scoef*dtit(i) + HII(i,j,k) )
     &                      / ( 1._DKIND + acoef*dtit(i) )
!     
!          3) electrons:

            scoef =   k8(i) * HM(i,j,k) * HI(i,j,k)
     &             +  k15(i)* HM(i,j,k) * HI(i,j,k)
     &             +  k17(i)* HM(i,j,k) * HII(i,j,k)
     &             +  k57(i)* HI(i,j,k) * HI(i,j,k)
     &             +  k58(i)* HI(i,j,k) * HeI(i,j,k)/4._DKIND
!                  
     &             + k24shield(i)*HIp(i)
     &             + k25shield(i)*HeIIp(i)/4._DKIND
     &             + k26shield(i)*HeIp(i)/4._DKIND

            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0) )
     &          scoef = scoef + kphHI(i,j,k) * HIp(i)
     &                + kphHeI(i,j,k)  * HeIp(i)  / 4._DKIND
     &                + kphHeII(i,j,k) * HeIIp(i) / 4._DKIND
            if ( (iradtrans .eq. 1) .and. (irt_honly .eq. 0) )
     &          scoef = scoef + kphHI(i,j,k) * HIp(i)

            acoef = - (k1(i) *HI(i,j,k)    - k2(i)*HII(i,j,k)
     &              +  k3(i) *HeI(i,j,k)/4._DKIND -
     &           k6(i)*HeIII(i,j,k)/4._DKIND
     &              +  k5(i) *HeII(i,j,k)/4._DKIND -
     &           k4(i)*HeII(i,j,k)/4._DKIND
     &              +  k14(i)*HM(i,j,k)
     &              -  k7(i) *HI(i,j,k)
     &              -  k18(i)*H2II(i,j,k)/2._DKIND)
            dep(i)  = ( scoef*dtit(i) + de(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

!           7) H2

            scoef = 2._DKIND*(k8(i)  * HM(i,j,k)   * HI(i,j,k)
     &            +       k10(i) * H2II(i,j,k) * HI(i,j,k)/2._DKIND
     &            +       k19(i) * H2II(i,j,k) * HM(i,j,k)/2._DKIND
     &            +       k22(i) * HI(i,j,k) * (HI(i,j,k))**2._DKIND)
            acoef = ( k13(i)*HI(i,j,k) + k11(i)*HII(i,j,k)
     &              + k12(i)*de(i,j,k) )
     &              + k29 + k31

            if (iradtrans .eq. 1) acoef = acoef + kdissH2I(i,j,k)

            if (idust .gt. 0) then
               scoef = scoef + 2._DKIND * h2dust(i) *
     &              HI(i,j,k) * rhoH(i)
            endif

            H2Ip(i) = ( scoef*dtit(i) + H2I(i,j,k) )
     &                / ( 1._DKIND + acoef*dtit(i) )

!           8) H-

            scoef = k7(i) * HI(i,j,k) * de(i,j,k) 
            acoef = (k8(i)  + k15(i))  * HI(i,j,k) + 
     &              (k16(i) + k17(i))  * HII(i,j,k) +  
     &	            k14(i) * de(i,j,k) + k19(i) * H2II(i,j,k)/2.0 +
     &	            k27
            HMp(i) = (scoef*dtit(i) + HM(i,j,k))
     &           / (1.0 + acoef*dtit(i))


!           9) H2+

            H2IIp(i) = 2._DKIND*( k9 (i)*HIp(i)*HIIp(i)
     &                    +   k11(i)*H2Ip(i)/2._DKIND*HIIp(i)
     &                    +   k17(i)*HMp(i)*HIIp(i)
     &                    + k29*H2Ip(i)
     &                    )
     &                 /  ( k10(i)*HIp(i) + k18(i)*dep(i)
     &                    + k19(i)*HMp(i)
     &                    + (k28+k30)
     &                    )

         endif                  ! itmask
         enddo
!     
      endif                     ! H2

!  --- (D) Now do extra 3-species for molecular HD ---
!     
      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if (itmask(i)) then
!     
!         1) DI
!     
            scoef =   (       k2(i) * DII(i,j,k) * de(i,j,k)
     &                 +      k51(i)* DII(i,j,k) * HI(i,j,k)
     &                 + 2._DKIND*k55(i)* HDI(i,j,k) *
     &              HI(i,j,k)/3._DKIND
     &                 )
            acoef  =    k1(i) * de(i,j,k)
     &             +    k50(i) * HII(i,j,k)
     &             +    k54(i) * H2I(i,j,k)/2._DKIND
     &             +    k56(i) * HM(i,j,k)
     &             + k24shield(i)
            if (iradtrans .eq. 1) acoef = acoef + kphHI(i,j,k)
            DIp(i)    = ( scoef*dtit(i) + DI(i,j,k) ) / 
     &                  ( 1._DKIND + acoef*dtit(i) )

!         2) DII
c 
            scoef =   (   k1(i)  * DI(i,j,k) * de(i,j,k)
     &            +       k50(i) * HII(i,j,k)* DI(i,j,k)
     &            +  2._DKIND*k53(i) * HII(i,j,k)* HDI(i,j,k)/3._DKIND
     &            )
     &            + k24shield(i)*DI(i,j,k)
            if (iradtrans .eq. 1) scoef = scoef + kphHI(i,j,k)*DI(i,j,k)
            acoef =    k2(i)  * de(i,j,k)
     &            +    k51(i) * HI(i,j,k)
     &            +    k52(i) * H2I(i,j,k)/2._DKIND

            DIIp(i)   = ( scoef*dtit(i) + DII(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

!          3) HDI
c 
            scoef = 3._DKIND*(k52(i) * DII(i,j,k)* 
     &           H2I(i,j,k)/2._DKIND/2._DKIND
     &           + k54(i) * DI(i,j,k) * H2I(i,j,k)/2._DKIND/2._DKIND
     &           + 2._DKIND*k56(i) * DI(i,j,k) * HM(i,j,k)/2._DKIND
     &                 )
            acoef  =    k53(i) * HII(i,j,k)
     &             +    k55(i) * HI(i,j,k)

            HDIp(i)   = ( scoef*dtit(i) + HDI(i,j,k) )
     &                 / ( 1._DKIND + acoef*dtit(i) )

         endif                  ! itmask
         enddo
      endif

!   --- (E) Set densities from 1D temps to 3D fields ---

      do i = is+1, ie+1
         if (itmask(i)) then
         HIdot_prev(i) = abs(HI(i,j,k)-HIp(i)) /
     &           max(real(dtit(i), DKIND), tiny8)
         HI(i,j,k)    = max(real(HIp(i), RKIND), tiny)
         HII(i,j,k)   = max(real(HIIp(i), RKIND), tiny)
         HeI(i,j,k)   = max(real(HeIp(i), RKIND), tiny)
         HeII(i,j,k)  = max(real(HeIIp(i), RKIND), tiny)
         HeIII(i,j,k) = max(real(HeIIIp(i), RKIND), 1e-5_RKIND*tiny)

!        de(i,j,k)    = dep(i)

!        Use charge conservation to determine electron fraction

         dedot_prev(i) = de(i,j,k)
         de(i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND +
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) 
     &        de(i,j,k) = de(i,j,k) - HM(i,j,k) + H2II(i,j,k)/2._RKIND
         dedot_prev(i) = abs(de(i,j,k)-dedot_prev(i))/
     &        max(dtit(i),tiny8)

         if (ispecies .gt. 1) then
            HM(i,j,k)    = max(real(HMp(i), RKIND), tiny)
            H2I(i,j,k)   = max(real(H2Ip(i),RKIND), tiny)
            H2II(i,j,k)  = max(real(H2IIp(i), RKIND), tiny)
         endif

         if (ispecies .gt. 2) then
            DI(i,j,k)    = max(real(DIp(i), RKIND), tiny)
            DII(i,j,k)   = max(real(DIIp(i), RKIND), tiny)
            HDI(i,j,k)   = max(real(HDIp(i), RKIND), tiny)
         endif
      endif                     ! itmask
!     

      if (HI(i,j,k) .ne. HI(i,j,k)) then
#ifdef _OPENMP
!$omp critical
         write(*,*) 'HUGE HI! :: ', i, j, k, HI(i,j,k)
!$omp end critical
#endif
      endif

      enddo                     ! end loop over i

      return
      end

! ------------------------------------------------------------------
!   This routine correct the highest abundence species to
!     insure conservation of particle number and charge.

      subroutine make_consistent_g(de, HI, HII, HeI, HeII, HeIII,
     &                        HM, H2I, H2II, DI, DII, HDI, metal, d,
     &                        is, ie, js, je, ks, ke,
     &                        in, jn, kn, ispecies, imetal, fh, dtoh)
! -------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"

!     Arguments

      integer in, jn, kn, is, ie, js, je, ks, ke, ispecies, imetal
      R_PREC  de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &        HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &        d(in,jn,kn), metal(in,jn,kn),
     &        HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn)
      real*8 fh, dtoh

!     locals

      integer i, j, k
      real*8 totalH(in), totalHe(in),
     &       totalD, metalfree(in)
      R_PREC correctH, correctHe, correctD

!     Loop over all zones

      do k = ks+1, ke+1
      do j = js+1, je+1

!     Compute total densities of H and He
!         (ensure non-negativity)

      if (imetal .eq. 1) then
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k) - metal(i,j,k)
         enddo
      else
         do i = is+1, ie+1
            metalfree(i) = d(i,j,k)
         enddo
      endif

      do i = is+1, ie+1
         HI   (i,j,k) = abs(HI   (i,j,k))
         HII  (i,j,k) = abs(HII  (i,j,k))
         HeI  (i,j,k) = abs(HeI  (i,j,k))
         HeII (i,j,k) = abs(HeII (i,j,k))
         HeIII(i,j,k) = abs(HeIII(i,j,k))
         totalH(i) = HI(i,j,k) + HII(i,j,k)
         totalHe(i) = HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k)
      enddo

!     include molecular hydrogen

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            HM   (i,j,k) = abs(HM   (i,j,k))
            H2II (i,j,k) = abs(H2II (i,j,k))
            H2I  (i,j,k) = abs(H2I  (i,j,k))
            totalH(i) = totalH(i) + HM(i,j,k) + H2I(i,j,k) + H2II(i,j,k)
         enddo
      endif

!     Correct densities by keeping fractions the same

      do i = is+1, ie+1
         correctH = real(fh*metalfree(i)/totalH(i), RKIND)
         HI(i,j,k)  = HI(i,j,k)*correctH
         HII(i,j,k) = HII(i,j,k)*correctH

         correctHe = real((1._DKIND - fh)*
     &        metalfree(i)/totalHe(i), RKIND)
         HeI(i,j,k)   = HeI(i,j,k)*correctHe
         HeII(i,j,k)  = HeII(i,j,k)*correctHe
         HeIII(i,j,k) = HeIII(i,j,k)*correctHe

!     Correct molecular hydrogen-related fractions

         if (ispecies .gt. 1) then
            HM   (i,j,k) = HM(i,j,k)*correctH
            H2II (i,j,k) = H2II(i,j,k)*correctH
            H2I  (i,j,k) = H2I(i,j,k)*correctH
         endif
      enddo

!     Do the same thing for deuterium (ignore HD) Assumes dtoh is small

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            DI  (i,j,k) = abs(DI  (i,j,k))
            DII (i,j,k) = abs(DII (i,j,k))
            HDI (i,j,k) = abs(HDI (i,j,k))
            totalD = DI(i,j,k) + DII(i,j,k) +
     &           2._DKIND/3._DKIND*HDI(i,j,k)
            correctD = real(fh*dtoh*metalfree(i)/totalD, RKIND)
            DI  (i,j,k) = DI (i,j,k)*correctD
            DII (i,j,k) = DII(i,j,k)*correctD
            HDI (i,j,k) = HDI(i,j,k)*correctD
         enddo
      endif

!       Set the electron density

      do i = is+1, ie+1
         de (i,j,k) = HII(i,j,k) + HeII(i,j,k)/4._RKIND +
     &        HeIII(i,j,k)/2._RKIND
         if (ispecies .gt. 1) de(i,j,k) = de(i,j,k)
     &        - HM(i,j,k) + H2II(i,j,k)/2._RKIND
      enddo

      enddo  ! end loop over j
      enddo  ! end loop over k

      return
      end
