#include "phys_const.def"

!=======================================================================
!\////////////////////  SUBROUTINE CALC_TDUST_3D  \\\\\\\\\\\\\\\\\\\\\\

      subroutine calc_tdust_3d_g(
     &     d, de, HI, HII, 
     &     HeI, HeII, HeIII,
     &     HM, H2I, H2II, 
     &     in, jn, kn, 
     &     nratec, iexpand,
     &     ispecies, idim,
     &     is, js, ks, 
     &     ie, je, ke, 
     &     aye, temstart, temend,
     &     fgr, gasgra,
     &     gamma_isrfa, isrf,
     &     utem, uxyz, uaye,
     &     urho, utim,
     &     gas_temp, dust_temp,
     &     iisrffield, isrf_habing
     &   , imetal, idustfield, igrgr
     &   , z_solar, metal, dust
     &   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &   , AC, SiO2D, MgO, FeS, Al2O3
     &   , reforg, volorg, H2Oice
     &   , immulti, impop3, idspecies, itdmulti, idsub
     &   , metal_loc, metal_C30, metal_F13
     &   , loc_fFeM, loc_fMg2SiO4, loc_fMgSiO3, loc_fFeS
     &   , loc_freforg, loc_fvolorg, loc_fH2Oice
     &   , loc_r0FeM, loc_r0Mg2SiO4, loc_r0MgSiO3, loc_r0FeS
     &   , loc_r0reforg, loc_r0volorg, loc_r0H2Oice 
     &   , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &   , C30_fAC, C30_fSiO2D, C30_fMgO
     &   , C30_fFeS, C30_fAl2O3
     &   , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4
     &   , C30_r0MgSiO3, C30_r0AC, C30_r0SiO2D
     &   , C30_r0MgO, C30_r0FeS, C30_r0Al2O3
     &   , F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &   , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fAl2O3
     &   , F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &   , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0Al2O3
     &   , gr_N, gr_Size, gr_dT, gr_Td
     &   , loc_kpFeM, loc_kpMg2SiO4, loc_kpMgSiO3, loc_kpFeS
     &   , loc_kpreforg, loc_kpvolorg, loc_kpH2Oice
     &   , C30_kpSiM, C30_kpFeM, C30_kpMg2SiO4
     &   , C30_kpMgSiO3, C30_kpAC, C30_kpSiO2D
     &   , C30_kpMgO, C30_kpFeS, C30_kpAl2O3
     &   , F13_kpFeM, F13_kpMg2SiO4, F13_kpMgSiO3
     &   , F13_kpFe3O4, F13_kpAC, F13_kpSiO2D, F13_kpAl2O3
     &   , SiM_temp, FeM_temp, Mg2SiO4_temp, MgSiO3_temp, Fe3O4_temp
     &   , AC_temp, SiO2D_temp, MgO_temp, FeS_temp, Al2O3_temp
     &   , reforg_temp, volorg_temp, H2Oice_temp
     &   , gasgr2a, gamma_isrf2a
     &   )

!  COMPUTE THE DUST TEMPERATURE
!
!  written by: Britton Smith
!  date: July, 2011
!  modified1: 
!
!  PURPOSE:
!    Calculate dust heat balance to get the dust temperature.
!
!  INPUTS:
!
!  PARAMETERS:
!
!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"
#ifdef _OPENMP
#include "omp_lib.h"
#endif

!  Arguments

      integer in, jn, kn, is, js, ks, ie, je, ke, nratec,
     &        iexpand, ispecies, idim, iisrffield
      real*8  aye, temstart, temend,
     &        utem, uxyz, uaye, urho, utim, fgr, isrf
      R_PREC  d(in,jn,kn),
     &     de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &     HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn),
     &     HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &     gas_temp(in,jn,kn), dust_temp(in,jn,kn),
     &     isrf_habing(in,jn,kn)
       real*8 dom
       real*8 z_solar, metallicity(in), dust2gas(in)
       R_PREC metal(in,jn,kn), dust(in,jn,kn)
       R_PREC metal_loc(in,jn,kn)
     &      , metal_C30(in,jn,kn), metal_F13(in,jn,kn)
     &   , CI(in,jn,kn)   , CII(in,jn,kn)    , CO(in,jn,kn)
     &   , CO2(in,jn,kn)  , OI(in,jn,kn)     , OH(in,jn,kn)
     &   , H2O(in,jn,kn)  , O2(in,jn,kn)     , SiI(in,jn,kn)
     &   , SiOI(in,jn,kn) , SiO2I(in,jn,kn)  , CH(in,jn,kn)
     &   , CH2(in,jn,kn)  , COII(in,jn,kn)   , OII(in,jn,kn)
     &   , OHII(in,jn,kn) , H2OII(in,jn,kn)  , H3OII(in,jn,kn)
     &   , O2II(in,jn,kn) , Mg(in,jn,kn)     , Al(in,jn,kn)
     &   , S(in,jn,kn)    , Fe(in,jn,kn)     
      R_PREC  SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn)
     &      , MgSiO3(in,jn,kn), Fe3O4(in,jn,kn), AC(in,jn,kn)
     &      , SiO2D(in,jn,kn), MgO(in,jn,kn), FeS(in,jn,kn)
     &      , Al2O3(in,jn,kn)
     &      , reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn)

!  Chemistry tables
      real*8  gasgra(nratec), gamma_isrfa
      integer imetal, idustfield, igrgr
      integer immulti, impop3, idspecies, itdmulti, idsub
      real*8  loc_fFeM, loc_fMg2SiO4, loc_fMgSiO3, loc_fFeS
     &      , loc_freforg, loc_fvolorg, loc_fH2Oice
      real*8  loc_r0FeM(3), loc_r0Mg2SiO4(3), loc_r0MgSiO3(3)
     &      , loc_r0FeS(3), loc_r0reforg(3), loc_r0volorg(3)
     &      , loc_r0H2Oice(3)
      real*8  C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &      , C30_fAC, C30_fSiO2D, C30_fMgO
     &      , C30_fFeS, C30_fAl2O3
      real*8  C30_r0SiM(3), C30_r0FeM(3), C30_r0Mg2SiO4(3)
     &      , C30_r0MgSiO3(3), C30_r0AC(3), C30_r0SiO2D(3)
     &      , C30_r0MgO(3), C30_r0FeS(3), C30_r0Al2O3(3)
      real*8  F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &      , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fAl2O3
      real*8  F13_r0FeM(3), F13_r0Mg2SiO4(3), F13_r0MgSiO3(3)
     &      , F13_r0Fe3O4(3), F13_r0AC(3), F13_r0SiO2D(3)
     &      , F13_r0Al2O3(3)
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  loc_kpFeM(gr_Size), loc_kpMg2SiO4(gr_Size)
     &      , loc_kpMgSiO3(gr_Size), loc_kpFeS(gr_Size)
     &      , loc_kpreforg(gr_Size), loc_kpvolorg(gr_Size)
     &      , loc_kpH2Oice(gr_Size)
      real*8  C30_kpSiM(gr_Size), C30_kpFeM(gr_Size)
     &      , C30_kpMg2SiO4(gr_Size), C30_kpMgSiO3(gr_Size)
     &      , C30_kpFe3O4(gr_Size), C30_kpAC(gr_Size)
     &      , C30_kpSiO2D(gr_Size), C30_kpMgO(gr_Size)
     &      , C30_kpFeS(gr_Size), C30_kpAl2O3(gr_Size)
      real*8  F13_kpFeM(gr_Size)
     &      , F13_kpMg2SiO4(gr_Size), F13_kpMgSiO3(gr_Size)
     &      , F13_kpFe3O4(gr_Size), F13_kpAC(gr_Size)
     &      , F13_kpSiO2D(gr_Size), F13_kpAl2O3(gr_Size)
!     grain growth
      real*8  sgSiM(in), sgFeM(in), sgMg2SiO4(in)
     &      , sgMgSiO3(in), sgFe3O4(in), sgAC(in)
     &      , sgSiO2D(in), sgMgO(in), sgFeS(in)
     &      , sgAl2O3(in)
     &      , sgreforg(in), sgvolorg(in), sgH2Oice(in)
     &      , sgtot(in)
      real*8  alSiM(gr_N(2),in), alFeM(gr_N(2),in)
     &      , alMg2SiO4(gr_N(2),in), alMgSiO3(gr_N(2),in)
     &      , alFe3O4(gr_N(2),in), alAC(gr_N(2),in)
     &      , alSiO2D(gr_N(2),in), alMgO(gr_N(2),in)
     &      , alFeS(gr_N(2),in), alAl2O3(gr_N(2),in)
     &      , alreforg(gr_N(2),in)
     &      , alvolorg(gr_N(2),in), alH2Oice(gr_N(2),in)
     &      , altot(gr_N(2),in)
      real*8  kpSiM(in), kpFeM(in), kpMg2SiO4(in)
     &      , kpMgSiO3(in), kpFe3O4(in), kpAC(in)
     &      , kpSiO2D(in), kpMgO(in), kpFeS(in)
     &      , kpAl2O3(in)
     &      , kpreforg(in), kpvolorg(in), kpH2Oice(in)
     &      , kptot(in)
!     grain temperature
      real*8  tSiM(in), tFeM(in), tMg2SiO4(in)
     &      , tMgSiO3(in), tFe3O4(in), tAC(in)
     &      , tSiO2D(in), tMgO(in), tFeS(in)
     &      , tAl2O3(in)
     &      , treforg(in), tvolorg(in), tH2Oice(in)
      R_PREC  SiM_temp(in,jn,kn), FeM_temp(in,jn,kn)
     &      , Mg2SiO4_temp(in,jn,kn), MgSiO3_temp(in,jn,kn)
     &      , Fe3O4_temp(in,jn,kn), AC_temp(in,jn,kn)
     &      , SiO2D_temp(in,jn,kn), MgO_temp(in,jn,kn)
     &      , FeS_temp(in,jn,kn), Al2O3_temp(in,jn,kn)
     &      , reforg_temp(in,jn,kn)
     &      , volorg_temp(in,jn,kn), H2Oice_temp(in,jn,kn)
      real*8  gasgr2a(nratec), gamma_isrf2a
      real*8  gasSiM(in), gasFeM(in), gasMg2SiO4(in)
     &      , gasMgSiO3(in), gasFe3O4(in), gasAC(in)
     &      , gasSiO2D(in), gasMgO(in), gasFeS(in)
     &      , gasAl2O3(in)
     &      , gasreforg(in), gasvolorg(in), gasH2Oice(in)
      real*8  gasgr_tSiM(in), gasgr_tFeM(in), gasgr_tMg2SiO4(in)
     &      , gasgr_tMgSiO3(in), gasgr_tFe3O4(in), gasgr_tAC(in)
     &      , gasgr_tSiO2D(in), gasgr_tMgO(in), gasgr_tFeS(in)
     &      , gasgr_tAl2O3(in)
     &      , gasgr_treforg(in), gasgr_tvolorg(in), gasgr_tH2Oice(in)
      real*8  mygisrf(in), fv2k, fac
      real*8  gisrfSiM(in), gisrfFeM(in), gisrfMg2SiO4(in)
     &      , gisrfMgSiO3(in), gisrfFe3O4(in), gisrfAC(in)
     &      , gisrfSiO2D(in), gisrfMgO(in), gisrfFeS(in)
     &      , gisrfAl2O3(in)
     &      , gisrfreforg(in), gisrfvolorg(in), gisrfH2Oice(in)

!  Parameters

      real*8 mh
      parameter (mh = mass_h)

!  Locals

      integer i, j, k
      integer t, dj, dk
      real*8  trad, zr, logtem0, logtem9, dlogtem,
     &        coolunit, dbase1, tbase1, xbase1

!  Slice locals
 
      integer indixe(in)
      R_PREC  t1(in), t2(in), logtem(in), tdef(in), 
     &     tgas(in), tdust(in), nh(in), gasgr(in),
     &     myisrf(in)
          
!  Iteration mask for multi_cool

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0  = log(temstart)
      logtem9  = log(temend)
      dlogtem  = (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      dom      = urho*(aye**3)/mh
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1._RKIND/(aye*uaye) - 1._RKIND

!     Set CMB temperature

      trad = 2.73_RKIND * (1._RKIND + zr)

!  Loop over zones, and do an entire i-column in one go
      dk = ke - ks + 1
      dj = je - js + 1

      if (iexpand .eq. 1) then

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   comp1, comp2, energy,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef, p2d,
!$omp&   tgas, tgasold,
!$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
!$omp&   mynh, myde, gammaha_eff, gasgr_tdust, regr, edot,
!$omp&   ceHI, ceHeI, ceHeII,
!$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
!$omp&   reHII, reHeII1, reHeII2, reHeIII,
!$omp&   brem, cieco,
!$omp&   hyd01k, h2k01, vibh, roth, rotl,
!$omp&   gpldl, gphdl, hdlte, hdlow,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
            if (idustfield .eq. 1) then
               do i = is+1, ie+1
                  dust(i,j,k) = dust(i,j,k)/real(aye**3, RKIND)
               enddo
            endif
            if (idspecies .eq. 1) then
!!!!!!!!!!     CHANGE-ME
            endif
      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif
 
      endif

!     Loop over slices (in the k-direction)

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   tgas, tdust, nh, gasgr, myisrf,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

!     Compute grain size increment

      if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then

         call calc_grain_size_increment_1d(
     &             igrgr, itmask
     &           , in, jn, kn, is, ie, j, k, dom, d
     &           , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &           , AC, SiO2D, MgO, FeS, Al2O3
     &           , reforg, volorg, H2Oice
     &           , metal, metal_loc, metal_C30, metal_F13
     &           , loc_fFeM, loc_fMg2SiO4, loc_fMgSiO3, loc_fFeS
     &           , loc_freforg, loc_fvolorg, loc_fH2Oice
     &           , loc_r0FeM, loc_r0Mg2SiO4, loc_r0MgSiO3, loc_r0FeS
     &           , loc_r0reforg, loc_r0volorg, loc_r0H2Oice 
     &           , C30_fSiM, C30_fFeM, C30_fMg2SiO4, C30_fMgSiO3
     &           , C30_fAC, C30_fSiO2D, C30_fMgO
     &           , C30_fFeS, C30_fAl2O3
     &           , C30_r0SiM, C30_r0FeM, C30_r0Mg2SiO4
     &           , C30_r0MgSiO3, C30_r0AC, C30_r0SiO2D
     &           , C30_r0MgO, C30_r0FeS, C30_r0Al2O3
     &           , F13_fFeM, F13_fMg2SiO4, F13_fMgSiO3
     &           , F13_fFe3O4, F13_fAC, F13_fSiO2D, F13_fAl2O3
     &           , F13_r0FeM, F13_r0Mg2SiO4, F13_r0MgSiO3
     &           , F13_r0Fe3O4, F13_r0AC, F13_r0SiO2D, F13_r0Al2O3
     &           , gr_N, gr_Size, gr_dT, gr_Td
     &           , loc_kpFeM, loc_kpMg2SiO4, loc_kpMgSiO3, loc_kpFeS
     &           , loc_kpreforg, loc_kpvolorg, loc_kpH2Oice
     &           , C30_kpSiM, C30_kpFeM, C30_kpMg2SiO4
     &           , C30_kpMgSiO3, C30_kpAC, C30_kpSiO2D
     &           , C30_kpMgO, C30_kpFeS, C30_kpAl2O3
     &           , F13_kpFeM, F13_kpMg2SiO4, F13_kpMgSiO3
     &           , F13_kpFe3O4, F13_kpAC, F13_kpSiO2D, F13_kpAl2O3
     &           , sgSiM, sgFeM, sgMg2SiO4, sgMgSiO3, sgFe3O4, sgAC
     &           , sgSiO2D, sgMgO, sgFeS, sgAl2O3
     &           , sgreforg, sgvolorg, sgH2Oice, sgtot
     &           , alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC
     &           , alSiO2D, alMgO, alFeS, alAl2O3
     &           , alreforg, alvolorg, alH2Oice, altot
     &         )

      endif

        do i = is+1, ie+1

!     Set itmask to all true

           itmask(i) = .true.

!     Compute interstellar radiation field

           if (iisrffield .gt. 0) then
              myisrf(i) = isrf_habing(i,j,k)
           else
              myisrf(i) = isrf
           endif

!     Calculate metallicity

           if (imetal .eq. 1) then
              metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
           endif

!     Calculate dust to gas ratio

           if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
               dust(i,j,k) =
     &             SiM     (i,j,k)
     &           + FeM     (i,j,k)
     &           + Mg2SiO4 (i,j,k)
     &           + MgSiO3  (i,j,k)
     &           + Fe3O4   (i,j,k)
     &           + AC      (i,j,k)
     &           + SiO2D   (i,j,k)
     &           + MgO     (i,j,k)
     &           + FeS     (i,j,k)
     &           + Al2O3   (i,j,k)
     &           + reforg  (i,j,k)
     &           + volorg  (i,j,k)
     &           + H2Oice  (i,j,k)
           endif

           if (idustfield .gt. 0) then
              dust2gas(i) = dust(i,j,k) / d(i,j,k)
           else
              dust2gas(i) = fgr * metallicity(i)
           endif

!     Compute interstellar radiation field

         if (idspecies .eq. 0 ) then
            if (idustfield .gt. 0) then
               mygisrf(i) = gamma_isrfa
     &                    * fgr / dust2gas(i) * metallicity(i)
               !! correct with the depletion or enhancement of condensation rate.
!!             write(*,*) 'a', mygisrf(i)
            else
               mygisrf(i) = gamma_isrfa
            endif

         else ! idspecies

            if (itdmulti .eq. 0) then
            
               mygisrf(i) = gamma_isrf2a * sgtot(i)
               !! in UV, absorption coefficient Q ~ 1 (Goldsmith 2001)
               !! so we use the geometrical cross-section of grains [cgs]
!!             write(*,*) 'b', mygisrf(i) / dust2gas(i)

            else

               gisrfSiM     (i) = gamma_isrf2a * sgSiM     (i)
               gisrfFeM     (i) = gamma_isrf2a * sgFeM     (i)
               gisrfMg2SiO4 (i) = gamma_isrf2a * sgMg2SiO4 (i)
               gisrfMgSiO3  (i) = gamma_isrf2a * sgMgSiO3  (i)
               gisrfFe3O4   (i) = gamma_isrf2a * sgFe3O4   (i)
               gisrfAC      (i) = gamma_isrf2a * sgAC      (i)
               gisrfSiO2D   (i) = gamma_isrf2a * sgSiO2D   (i)
               gisrfMgO     (i) = gamma_isrf2a * sgMgO     (i)
               gisrfFeS     (i) = gamma_isrf2a * sgFeS     (i)
               gisrfAl2O3   (i) = gamma_isrf2a * sgAl2O3   (i)
               gisrfreforg  (i) = gamma_isrf2a * sgreforg  (i)
               gisrfvolorg  (i) = gamma_isrf2a * sgvolorg  (i)
               gisrfH2Oice  (i) = gamma_isrf2a * sgH2Oice  (i)

            endif

         endif

!     Compute hydrogen number density

           nh(i) = HI(i,j,k) + HII(i,j,k)
           if (ispecies .gt. 1) then
              nh(i) = nh(i) + H2I(i,j,k) + H2II(i,j,k)
           endif

!     We have not converted to proper, so use urho and not dom

           nh(i) = nh(i) * urho / mh

!     Compute log temperature and truncate if above/below table max/min

           tgas(i)   = gas_temp(i,j,k)
           logtem(i) = log(tgas(i))
           logtem(i) = max(logtem(i), logtem0)
           logtem(i) = min(logtem(i), logtem9)

!     Compute index into the table and precompute parts of linear interp

           indixe(i) = min(nratec-1,
     &          max(1,int((logtem(i)-logtem0)/dlogtem,DIKIND)+1))
           t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
           t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
           tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!     Lookup values and do a linear temperature in log(T)
!     Convert back to cgs

            if(idspecies .eq. 0) then

               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))

!!             gasgr(i) = fgr * gasgr(i) * coolunit / mh
               gasgr(i) = (dust2gas(i) / metallicity(i))
     &                              * gasgr(i) * coolunit / mh
                    !! apply to (idustfield .eq. 1) GC20200701
!!             write(*,*) 'a', gasgr(i)

            else

               fv2k    = gasgr2a(indixe(i)) + tdef(i)
     &              *(gasgr2a(indixe(i)+1) -gasgr2a(indixe(i)))
                 
               fac = coolunit / mh

               if ( itdmulti .eq. 0 ) then

               gasgr(i) = fv2k * sgtot(i)

               gasgr(i) = gasgr(i) * dust(i,j,k) * fac
!!             write(*,*) 'b', gasgr(i) / metallicity(i)

                 
               else

               gasSiM     (i) = fv2k * sgSiM     (i)
               gasFeM     (i) = fv2k * sgFeM     (i)
               gasMg2SiO4 (i) = fv2k * sgMg2SiO4 (i)
               gasMgSiO3  (i) = fv2k * sgMgSiO3  (i)
               gasFe3O4   (i) = fv2k * sgFe3O4   (i)
               gasAC      (i) = fv2k * sgAC      (i)
               gasSiO2D   (i) = fv2k * sgSiO2D   (i)
               gasMgO     (i) = fv2k * sgMgO     (i)
               gasFeS     (i) = fv2k * sgFeS     (i)
               gasAl2O3   (i) = fv2k * sgAl2O3   (i)
               gasreforg  (i) = fv2k * sgreforg  (i)
               gasvolorg  (i) = fv2k * sgvolorg  (i)
               gasH2Oice  (i) = fv2k * sgH2Oice  (i)

               gasSiM     (i) = gasSiM     (i) * fac
               gasFeM     (i) = gasFeM     (i) * fac
               gasMg2SiO4 (i) = gasMg2SiO4 (i) * fac
               gasMgSiO3  (i) = gasMgSiO3  (i) * fac
               gasFe3O4   (i) = gasFe3O4   (i) * fac
               gasAC      (i) = gasAC      (i) * fac
               gasSiO2D   (i) = gasSiO2D   (i) * fac
               gasMgO     (i) = gasMgO     (i) * fac
               gasFeS     (i) = gasFeS     (i) * fac
               gasAl2O3   (i) = gasAl2O3   (i) * fac
               gasreforg  (i) = gasreforg  (i) * fac
               gasvolorg  (i) = gasvolorg  (i) * fac
               gasH2Oice  (i) = gasH2Oice  (i) * fac

               endif

            endif


        enddo

!     --- Compute dust temperature in a slice ---

         if (itdmulti .eq. 0) then

            call calc_tdust_1d_g(tdust, tgas, nh, gasgr,
     &           mygisrf, myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, altot, kptot)

!!          write(*,*) tdust

         else

            call calc_tdust_1d_g(tSiM     , tgas, nh, gasgr,
     &           gisrfSiM     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiM     , kpSiM     )
          
            call calc_tdust_1d_g(tFeM     , tgas, nh, gasgr,
     &           gisrfFeM     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeM     , kpFeM     )
          
            call calc_tdust_1d_g(tMg2SiO4 , tgas, nh, gasgr,
     &           gisrfMg2SiO4 , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMg2SiO4 , kpMg2SiO4 )
          
            call calc_tdust_1d_g(tMgSiO3  , tgas, nh, gasgr,
     &           gisrfMgSiO3  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgSiO3  , kpMgSiO3  )
          
            call calc_tdust_1d_g(tFe3O4   , tgas, nh, gasgr,
     &           gisrfFe3O4   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFe3O4   , kpFe3O4   )
          
            call calc_tdust_1d_g(tAC      , tgas, nh, gasgr,
     &           gisrfAC      , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAC      , kpAC      )
          
            call calc_tdust_1d_g(tSiO2D   , tgas, nh, gasgr,
     &           gisrfSiO2D   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiO2D   , kpSiO2D   )
          
            call calc_tdust_1d_g(tMgO     , tgas, nh, gasgr,
     &           gisrfMgO     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgO     , kpMgO     )
          
            call calc_tdust_1d_g(tFeS     , tgas, nh, gasgr,
     &           gisrfFeS     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeS     , kpFeS     )
          
            call calc_tdust_1d_g(tAl2O3   , tgas, nh, gasgr,
     &           gisrfAl2O3   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAl2O3   , kpAl2O3   )
          
            call calc_tdust_1d_g(treforg  , tgas, nh, gasgr,
     &           gisrfreforg  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alreforg  , kpreforg  )
          
            call calc_tdust_1d_g(tvolorg  , tgas, nh, gasgr,
     &           gisrfvolorg  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alvolorg  , kpvolorg  )
          
            call calc_tdust_1d_g(tH2Oice  , tgas, nh, gasgr,
     &           gisrfH2Oice  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alH2Oice  , kpH2Oice  )

         endif

!     Copy slice values back to grid

        do i = is+1, ie+1
           if (itdmulti .gt. 0) then
                SiM_temp(i,j,k) = tSiM     (i)
                FeM_temp(i,j,k) = tFeM     (i)
            Mg2SiO4_temp(i,j,k) = tMg2SiO4 (i)
             MgSiO3_temp(i,j,k) = tMgSiO3  (i)
              Fe3O4_temp(i,j,k) = tFe3O4   (i)
                 AC_temp(i,j,k) = tAC      (i)
              SiO2D_temp(i,j,k) = tSiO2D   (i)
                MgO_temp(i,j,k) = tMgO     (i)
                FeS_temp(i,j,k) = tFeS     (i)
              Al2O3_temp(i,j,k) = tAl2O3   (i)
             reforg_temp(i,j,k) = treforg  (i)
             volorg_temp(i,j,k) = tvolorg  (i)
             H2Oice_temp(i,j,k) = tH2Oice  (i)
           else
              dust_temp(i,j,k) = tdust(i)
           endif
        enddo

      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

!     Convert densities back to comoving from 'proper'

      if (iexpand .eq. 1) then

! parallelize the k and j loops with OpenMP
! flat j and k loops for better parallelism
#ifdef _OPENMP
!$omp parallel do schedule(runtime) private(
!$omp&   i, j, k,
!$omp&   comp1, comp2, energy,
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef, p2d,
!$omp&   tgas, tgasold,
!$omp&   tdust, metallicity, dust2gas, rhoH, mmw,
!$omp&   mynh, myde, gammaha_eff, gasgr_tdust, regr, edot,
!$omp&   ceHI, ceHeI, ceHeII,
!$omp&   ciHI, ciHeI, ciHeIS, ciHeII,
!$omp&   reHII, reHeII1, reHeII2, reHeIII,
!$omp&   brem, cieco,
!$omp&   hyd01k, h2k01, vibh, roth, rotl,
!$omp&   gpldl, gphdl, hdlte, hdlow,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

            if (imetal .eq. 1) then
               do i = is+1, ie+1
                  metal(i,j,k) = metal(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
            if (idustfield .eq. 1) then
               do i = is+1, ie+1
                  dust(i,j,k) = dust(i,j,k)*real(aye**3, RKIND)
               enddo
            endif
            if (idspecies .eq. 1) then
!!!!!!!!!      CHANGE-ME
            endif
      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

      endif

      return
      end
