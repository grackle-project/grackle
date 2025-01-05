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
     &   , imetal, imchem, idustfield, igrgr
     &   , z_solar, metal, dust
     &   , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &   , AC, SiO2D, MgO, FeS, Al2O3
     &   , reforg, volorg, H2Oice
     &   , immulti, imabund, idspecies, itdmulti, idsub
     &   , metal_loc
     &   , metal_C13, metal_C20, metal_C25, metal_C30
     &   , metal_F13, metal_F15, metal_F50, metal_F80
     &   , metal_P170, metal_P200, metal_Y19
     &   , SN0_N
     &   , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &   , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &   , SN0_fFeS, SN0_fAl2O3
     &   , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &   , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &   , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &   , SN0_r0FeS, SN0_r0Al2O3
     &   , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &   , gr_N, gr_Size, gr_dT, gr_Td
     &   , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &   , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &   , SN0_kpFeS, SN0_kpAl2O3
     &   , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &   , gasgr2a, gamma_isrf2a
     &   , SiM_temp, FeM_temp, Mg2SiO4_temp, MgSiO3_temp, Fe3O4_temp
     &   , AC_temp, SiO2D_temp, MgO_temp, FeS_temp, Al2O3_temp
     &   , reforg_temp, volorg_temp, H2Oice_temp
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
     &      , metal_C13(in,jn,kn), metal_C20(in,jn,kn)
     &      , metal_C25(in,jn,kn), metal_C30(in,jn,kn)
     &      , metal_F13(in,jn,kn), metal_F15(in,jn,kn)
     &      , metal_F50(in,jn,kn), metal_F80(in,jn,kn)
     &      , metal_P170(in,jn,kn), metal_P200(in,jn,kn)
     &      , metal_Y19(in,jn,kn)
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
      integer imetal, idustfield, imchem, igrgr
      integer immulti, imabund, idspecies, itdmulti, idsub
      integer SN0_N
      real*8  SN0_fSiM(SN0_N), SN0_fFeM(SN0_N), SN0_fMg2SiO4(SN0_N)
     &      , SN0_fMgSiO3(SN0_N), SN0_fFe3O4(SN0_N), SN0_fAC(SN0_N)
     &      , SN0_fSiO2D(SN0_N), SN0_fMgO(SN0_N), SN0_fFeS(SN0_N)
     &      , SN0_fAl2O3(SN0_N)
     &      , SN0_freforg(SN0_N), SN0_fvolorg(SN0_N), SN0_fH2Oice(SN0_N)
      real*8  SN0_r0SiM(3,SN0_N), SN0_r0FeM(3,SN0_N)
     &      , SN0_r0Mg2SiO4(3,SN0_N), SN0_r0MgSiO3(3,SN0_N)
     &      , SN0_r0Fe3O4(3,SN0_N), SN0_r0AC(3,SN0_N)
     &      , SN0_r0SiO2D(3,SN0_N), SN0_r0MgO(3,SN0_N)
     &      , SN0_r0FeS(3,SN0_N), SN0_r0Al2O3(3,SN0_N)
     &      , SN0_r0reforg(3,SN0_N)
     &      , SN0_r0volorg(3,SN0_N), SN0_r0H2Oice(3,SN0_N)
!     opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN0_kpSiM(gr_Size,SN0_N), SN0_kpFeM(gr_Size,SN0_N)
     &      , SN0_kpMg2SiO4(gr_Size,SN0_N), SN0_kpMgSiO3(gr_Size,SN0_N)
     &      , SN0_kpFe3O4(gr_Size,SN0_N), SN0_kpAC(gr_Size,SN0_N)
     &      , SN0_kpSiO2D(gr_Size,SN0_N), SN0_kpMgO(gr_Size,SN0_N)
     &      , SN0_kpFeS(gr_Size,SN0_N), SN0_kpAl2O3(gr_Size,SN0_N)
     &      , SN0_kpreforg(gr_Size,SN0_N)
     &      , SN0_kpvolorg(gr_Size,SN0_N), SN0_kpH2Oice(gr_Size,SN0_N)
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

      real*8, parameter :: mh = mass_h

!  Locals

      integer i, j, k
      integer t, dj, dk
      real*8  trad, zr, logtem0, logtem9, dlogtem,
     &        coolunit, dbase1, tbase1, xbase1
      integer itd

!  Slice locals
 
      integer indixe(in)
      R_PREC  t1(in), t2(in), logtem(in), tdef(in), 
     &     tgas(in), tdust(in), nh(in), gasgr(in),
     &     myisrf(in)

!  Iteration mask for multi_cool

      MASK_TYPE itmask(in)

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
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   tgas,
!$omp&   tdust, metallicity, dust2gas,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

        if (imetal .eq. 1) then
           do i = is+1, ie+1
              metal(i,j,k) = metal(i,j,k)/real(aye**3, RKIND)
              if (immulti .gt. 0) then
                 metal_loc(i,j,k) = metal_loc(i,j,k)/real(aye**3, RKIND)
                 metal_C13(i,j,k) = metal_C13(i,j,k)/real(aye**3, RKIND)
                 metal_C20(i,j,k) = metal_C20(i,j,k)/real(aye**3, RKIND)
                 metal_C25(i,j,k) = metal_C25(i,j,k)/real(aye**3, RKIND)
                 metal_C30(i,j,k) = metal_C30(i,j,k)/real(aye**3, RKIND)
                 metal_F13(i,j,k) = metal_F13(i,j,k)/real(aye**3, RKIND)
                 metal_F15(i,j,k) = metal_F15(i,j,k)/real(aye**3, RKIND)
                 metal_F50(i,j,k) = metal_F50(i,j,k)/real(aye**3, RKIND)
                 metal_F80(i,j,k) = metal_F80(i,j,k)/real(aye**3, RKIND)
                 metal_P170(i,j,k)=metal_P170(i,j,k)/real(aye**3, RKIND)
                 metal_P200(i,j,k)=metal_P200(i,j,k)/real(aye**3, RKIND)
                 metal_Y19(i,j,k) = metal_Y19(i,j,k)/real(aye**3, RKIND)
              endif
!!            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (imchem .eq. 1) then
                 CI(i,j,k)      = CI(i,j,k)/real(aye**3, RKIND)
                 CII(i,j,k)     = CII(i,j,k)/real(aye**3, RKIND)
                 CO(i,j,k)      = CO(i,j,k)/real(aye**3, RKIND)
                 CO2(i,j,k)     = CO2(i,j,k)/real(aye**3, RKIND)
                 OI(i,j,k)      = OI(i,j,k)/real(aye**3, RKIND)
                 OH(i,j,k)      = OH(i,j,k)/real(aye**3, RKIND)
                 H2O(i,j,k)     = H2O(i,j,k)/real(aye**3, RKIND)
                 O2(i,j,k)      = O2(i,j,k)/real(aye**3, RKIND)
                 SiI(i,j,k)     = SiI(i,j,k)/real(aye**3, RKIND)
                 SiOI(i,j,k)    = SiOI(i,j,k)/real(aye**3, RKIND)
                 SiO2I(i,j,k)   = SiO2I(i,j,k)/real(aye**3, RKIND)
                 CH(i,j,k)      = CH(i,j,k)/real(aye**3, RKIND)
                 CH2(i,j,k)     = CH2(i,j,k)/real(aye**3, RKIND)
                 COII(i,j,k)    = COII(i,j,k)/real(aye**3, RKIND)
                 OII(i,j,k)     = OII(i,j,k)/real(aye**3, RKIND)
                 OHII(i,j,k)    = OHII(i,j,k)/real(aye**3, RKIND)
                 H2OII(i,j,k)   = H2OII(i,j,k)/real(aye**3, RKIND)
                 H3OII(i,j,k)   = H3OII(i,j,k)/real(aye**3, RKIND)
                 O2II(i,j,k)    = O2II(i,j,k)/real(aye**3, RKIND)
                 if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                 if (idspecies .gt. 0) then
                 Mg(i,j,k)      = Mg(i,j,k)/real(aye**3, RKIND)
                 endif
                 if (idspecies .gt. 1) then
                 Al(i,j,k)      = Al(i,j,k)/real(aye**3, RKIND)
                 S(i,j,k)       = S(i,j,k)/real(aye**3, RKIND)
                 Fe(i,j,k)      = Fe(i,j,k)/real(aye**3, RKIND)
                 endif
                 endif
              endif
!!            endif
           enddo
        endif
        if (idustfield .eq. 1) then
           do i = is+1, ie+1
              dust(i,j,k) = dust(i,j,k)/real(aye**3, RKIND)
              if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
!!            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (idspecies .gt. 0) then
                 MgSiO3(i,j,k)  = MgSiO3(i,j,k)/real(aye**3, RKIND)
                 AC(i,j,k)      = AC(i,j,k)/real(aye**3, RKIND)
              endif
              if (idspecies .gt. 1) then
                 SiM(i,j,k)     = SiM(i,j,k)/real(aye**3, RKIND)
                 FeM(i,j,k)     = FeM(i,j,k)/real(aye**3, RKIND)
                 Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)/real(aye**3, RKIND)
                 Fe3O4(i,j,k)   = Fe3O4(i,j,k)/real(aye**3, RKIND)
                 SiO2D(i,j,k)   = SiO2D(i,j,k)/real(aye**3, RKIND)
                 MgO(i,j,k)     = MgO(i,j,k)/real(aye**3, RKIND)
                 FeS(i,j,k)     = FeS(i,j,k)/real(aye**3, RKIND)
                 Al2O3(i,j,k)   = Al2O3(i,j,k)/real(aye**3, RKIND)
              endif
              if (idspecies .gt. 2) then
                 reforg(i,j,k)  = reforg(i,j,k)/real(aye**3, RKIND)
                 volorg(i,j,k)  = volorg(i,j,k)/real(aye**3, RKIND)
                 H2Oice(i,j,k)  = H2Oice(i,j,k)/real(aye**3, RKIND)
              endif
!!            endif
              endif
           enddo
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

        do i = is+1, ie+1

!     Set itmask to all true

           itmask(i) = MASK_TRUE

        enddo

!     Iteration mask for metal-rich cells

         if (imetal .eq. 1) then
             do i = is+1, ie + 1
                if (metal(i,j,k) .lt. 1.e-9_DKIND * d(i,j,k)) then
                   itmask(i) = MASK_FALSE
                endif
             enddo
         endif

!     Compute grain size increment

      if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then

         call calc_grain_size_increment_1d(
     &             immulti, imabund, idspecies, igrgr, itmask
     &           , in, jn, kn, is, ie, j, k, dom, d
     &           , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &           , AC, SiO2D, MgO, FeS, Al2O3
     &           , reforg, volorg, H2Oice
     &           , metal, metal_loc
     &           , metal_C13, metal_C20, metal_C25, metal_C30
     &           , metal_F13, metal_F15, metal_F50, metal_F80
     &           , metal_P170, metal_P200, metal_Y19
     &           , SN0_N
     &           , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &           , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &           , SN0_fFeS, SN0_fAl2O3
     &           , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &           , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &           , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &           , SN0_r0FeS, SN0_r0Al2O3
     &           , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &           , gr_N, gr_Size, gr_dT, gr_Td
     &           , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &           , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &           , SN0_kpFeS, SN0_kpAl2O3
     &           , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &           , sgSiM, sgFeM, sgMg2SiO4, sgMgSiO3, sgFe3O4, sgAC
     &           , sgSiO2D, sgMgO, sgFeS, sgAl2O3
     &           , sgreforg, sgvolorg, sgH2Oice, sgtot
     &           , alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC
     &           , alSiO2D, alMgO, alFeS, alAl2O3
     &           , alreforg, alvolorg, alH2Oice, altot
     &         )

        endif

        do i = is+1, ie+1
        if(itmask(i) .ne. MASK_FALSE) then
!     Calculate metallicity

           if (imetal .eq. 1) then
              metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
           endif

!     Calculate dust to gas ratio

C          if ( (idustfield .gt. 0) .and. (idspecies .gt. 0) ) then
C           if (idspecies .gt. 0) then
C              dust(i,j,k) = MgSiO3  (i,j,k)
C    &                     + AC      (i,j,k)
C           endif
C           if (idspecies .gt. 1) then
C              dust(i,j,k) = dust(i,j,k)
C    &                     + SiM     (i,j,k)
C    &                     + FeM     (i,j,k)
C    &                     + Mg2SiO4 (i,j,k)
C    &                     + Fe3O4   (i,j,k)
C    &                     + SiO2D   (i,j,k)
C    &                     + MgO     (i,j,k)
C    &                     + FeS     (i,j,k)
C    &                     + Al2O3   (i,j,k)
C           endif
C           if (idspecies .gt. 2) then
C              dust(i,j,k) = dust(i,j,k)
C    &                     + reforg  (i,j,k)
C    &                     + volorg  (i,j,k)
C    &                     + H2Oice  (i,j,k)
C           endif
C          endif

           if (idustfield .gt. 0) then
              dust2gas(i) = dust(i,j,k) / d(i,j,k)
           else
              dust2gas(i) = fgr * metallicity(i)
           endif

!     Compute interstellar radiation field

           if (iisrffield .gt. 0) then
              myisrf(i) = isrf_habing(i,j,k)
           else
              myisrf(i) = isrf
           endif

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

            if (idspecies .gt. 0) then
               gisrfMgSiO3  (i) = gamma_isrf2a * sgMgSiO3  (i)
!!             write(*,*) d(i,j,k), gamma_isrf2a, sgMgSiO3(i)
               gisrfAC      (i) = gamma_isrf2a * sgAC      (i)
            endif
            if (idspecies .gt. 1) then
               gisrfSiM     (i) = gamma_isrf2a * sgSiM     (i)
               gisrfFeM     (i) = gamma_isrf2a * sgFeM     (i)
               gisrfMg2SiO4 (i) = gamma_isrf2a * sgMg2SiO4 (i)
               gisrfFe3O4   (i) = gamma_isrf2a * sgFe3O4   (i)
               gisrfSiO2D   (i) = gamma_isrf2a * sgSiO2D   (i)
               gisrfMgO     (i) = gamma_isrf2a * sgMgO     (i)
               gisrfFeS     (i) = gamma_isrf2a * sgFeS     (i)
               gisrfAl2O3   (i) = gamma_isrf2a * sgAl2O3   (i)
            endif
            if (idspecies .gt. 2) then
               gisrfreforg  (i) = gamma_isrf2a * sgreforg  (i)
               gisrfvolorg  (i) = gamma_isrf2a * sgvolorg  (i)
               gisrfH2Oice  (i) = gamma_isrf2a * sgH2Oice  (i)
            endif

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

               gasgr(i) = fv2k * sgtot(i) * fac
!!             write(*,*) 'b', gasgr(i) / metallicity(i)

               else

               if (idspecies .gt. 0) then
               gasMgSiO3  (i) = fv2k * sgMgSiO3  (i) * fac
!!             write(*,*) fv2k, sgMgSiO3(i), fac
               gasAC      (i) = fv2k * sgAC      (i) * fac
               endif
               if (idspecies .gt. 1) then
               gasSiM     (i) = fv2k * sgSiM     (i) * fac
               gasFeM     (i) = fv2k * sgFeM     (i) * fac
               gasMg2SiO4 (i) = fv2k * sgMg2SiO4 (i) * fac
               gasFe3O4   (i) = fv2k * sgFe3O4   (i) * fac
               gasSiO2D   (i) = fv2k * sgSiO2D   (i) * fac
               gasMgO     (i) = fv2k * sgMgO     (i) * fac
               gasFeS     (i) = fv2k * sgFeS     (i) * fac
               gasAl2O3   (i) = fv2k * sgAl2O3   (i) * fac
               endif
               if (idspecies .gt. 2) then
               gasreforg  (i) = fv2k * sgreforg  (i) * fac
               gasvolorg  (i) = fv2k * sgvolorg  (i) * fac
               gasH2Oice  (i) = fv2k * sgH2Oice  (i) * fac
               endif

               endif

            endif

        endif ! itmask
        enddo

!        if (idspecies .eq. 0) then
!           do itd = 1, gr_N(2)
!              do i = is+1, ie + 1
!              if(itmask(i)) then
!                 altot(itd,i)=
!    &                (SN0_kpFeM    (4*itd,1)*SN0_fFeM    (1)
!    &                /SN0_r0FeM    (3    ,1)
!    &                +SN0_kpMg2SiO4(4*itd,1)*SN0_fMg2SiO4(1)
!    &                /SN0_r0Mg2SiO4(3    ,1)
!    &                +SN0_kpMgSiO3 (4*itd,1)*SN0_fMgSiO3 (1)
!    &                /SN0_r0MgSiO3 (3    ,1)
!    &                +SN0_kpFeS    (4*itd,1)*SN0_fFeS    (1)
!    &                /SN0_r0FeS    (3    ,1)
!    &                +SN0_kpreforg (4*itd,1)*SN0_freforg (1)
!    &                /SN0_r0reforg (3    ,1)
!    &                +SN0_kpvolorg (4*itd,1)*SN0_fvolorg (1)
!    &                /SN0_r0volorg (3    ,1)
!    &                +SN0_kpH2Oice (4*itd,1)*SN0_fH2Oice (1)
!    &                /SN0_r0H2Oice (3    ,1))
!    &                  * metal(i,j,k) / d(i,j,k)
!              endif
!              enddo
!           enddo
!        endif

!     --- Compute dust temperature in a slice ---

         if (itdmulti .eq. 0) then

            call calc_tdust_1d_g(tdust, tgas, nh, gasgr,
     &           mygisrf, myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, altot, kptot, idspecies)

!!          write(*,*) tdust

         else
          
            if (idspecies .gt. 0) then
            call calc_tdust_1d_g(tMgSiO3  , tgas, nh, gasMgSiO3  ,
     &           gisrfMgSiO3  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgSiO3  , kpMgSiO3 
     &         , idspecies)
!!          write(*,*) 'bb2', gasMgSiO3, gisrfMgSiO3, myisrf, tMgSiO3
          
            call calc_tdust_1d_g(tAC      , tgas, nh, gasAC      ,
     &           gisrfAC      , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAC      , kpAC       
     &         , idspecies)
!!          write(*,*) 'bb2', gasAC    , gisrfAC    , myisrf, tAC
            endif

            if (idspecies .gt. 1) then
            call calc_tdust_1d_g(tSiM     , tgas, nh, gasSiM     ,
     &           gisrfSiM     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiM     , kpSiM      
     &         , idspecies)
          
            call calc_tdust_1d_g(tFeM     , tgas, nh, gasFeM     ,
     &           gisrfFeM     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeM     , kpFeM      
     &         , idspecies)
          
            call calc_tdust_1d_g(tMg2SiO4 , tgas, nh, gasMg2SiO4 ,
     &           gisrfMg2SiO4 , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMg2SiO4 , kpMg2SiO4  
     &         , idspecies)
          
            call calc_tdust_1d_g(tFe3O4   , tgas, nh, gasFe3O4   ,
     &           gisrfFe3O4   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFe3O4   , kpFe3O4    
     &         , idspecies)
          
            call calc_tdust_1d_g(tSiO2D   , tgas, nh, gasSiO2D   ,
     &           gisrfSiO2D   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiO2D   , kpSiO2D    
     &         , idspecies)
          
            call calc_tdust_1d_g(tMgO     , tgas, nh, gasMgO     ,
     &           gisrfMgO     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgO     , kpMgO      
     &         , idspecies)
          
            call calc_tdust_1d_g(tFeS     , tgas, nh, gasFeS     ,
     &           gisrfFeS     , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeS     , kpFeS      
     &         , idspecies)
          
            call calc_tdust_1d_g(tAl2O3   , tgas, nh, gasAl2O3   ,
     &           gisrfAl2O3   , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAl2O3   , kpAl2O3    
     &         , idspecies)
            endif

            if (idspecies .gt. 2) then
            call calc_tdust_1d_g(treforg  , tgas, nh, gasreforg  ,
     &           gisrfreforg  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alreforg  , kpreforg   
     &         , idspecies)
          
            call calc_tdust_1d_g(tvolorg  , tgas, nh, gasvolorg  ,
     &           gisrfvolorg  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alvolorg  , kpvolorg   
     &         , idspecies)
          
            call calc_tdust_1d_g(tH2Oice  , tgas, nh, gasH2Oice  ,
     &           gisrfH2Oice  , myisrf, itmask, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alH2Oice  , kpH2Oice   
     &         , idspecies)
            endif

         endif

!     Copy slice values back to grid

        do i = is+1, ie+1
        if (itmask(i) .ne. MASK_FALSE) then
           if (itdmulti .eq. 0) then
              dust_temp(i,j,k) = tdust(i)
           else
             if (idspecies .gt. 0) then
             MgSiO3_temp(i,j,k) = tMgSiO3  (i)
                 AC_temp(i,j,k) = tAC      (i)
             endif
             if (idspecies .gt. 1) then
                SiM_temp(i,j,k) = tSiM     (i)
                FeM_temp(i,j,k) = tFeM     (i)
            Mg2SiO4_temp(i,j,k) = tMg2SiO4 (i)
              Fe3O4_temp(i,j,k) = tFe3O4   (i)
              SiO2D_temp(i,j,k) = tSiO2D   (i)
                MgO_temp(i,j,k) = tMgO     (i)
                FeS_temp(i,j,k) = tFeS     (i)
              Al2O3_temp(i,j,k) = tAl2O3   (i)
             endif
             if (idspecies .gt. 2) then
             reforg_temp(i,j,k) = treforg  (i)
             volorg_temp(i,j,k) = tvolorg  (i)
             H2Oice_temp(i,j,k) = tH2Oice  (i)
             endif
           endif
        endif ! itmask
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
!$omp&   indixe,
!$omp&   t1, t2, logtem, tdef,
!$omp&   tgas,
!$omp&   tdust, metallicity, dust2gas,
!$omp&   itmask )
#endif
      do t = 0, dk*dj-1
        k = t/dj      + ks+1
        j = mod(t,dj) + js+1

        if (imetal .eq. 1) then
           do i = is+1, ie+1
              metal(i,j,k) = metal(i,j,k)*real(aye**3, RKIND)
              if (immulti .gt. 0) then
                 metal_loc(i,j,k) = metal_loc(i,j,k)*real(aye**3, RKIND)
                 metal_C13(i,j,k) = metal_C13(i,j,k)*real(aye**3, RKIND)
                 metal_C20(i,j,k) = metal_C20(i,j,k)*real(aye**3, RKIND)
                 metal_C25(i,j,k) = metal_C25(i,j,k)*real(aye**3, RKIND)
                 metal_C30(i,j,k) = metal_C30(i,j,k)*real(aye**3, RKIND)
                 metal_F13(i,j,k) = metal_F13(i,j,k)*real(aye**3, RKIND)
                 metal_F15(i,j,k) = metal_F15(i,j,k)*real(aye**3, RKIND)
                 metal_F50(i,j,k) = metal_F50(i,j,k)*real(aye**3, RKIND)
                 metal_F80(i,j,k) = metal_F80(i,j,k)*real(aye**3, RKIND)
                 metal_P170(i,j,k)=metal_P170(i,j,k)*real(aye**3, RKIND)
                 metal_P200(i,j,k)=metal_P200(i,j,k)*real(aye**3, RKIND)
                 metal_Y19(i,j,k) = metal_Y19(i,j,k)*real(aye**3, RKIND)
              endif
!!            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (imchem .eq. 1) then
                 CI(i,j,k)      = CI(i,j,k)*real(aye**3, RKIND)
                 CII(i,j,k)     = CII(i,j,k)*real(aye**3, RKIND)
                 CO(i,j,k)      = CO(i,j,k)*real(aye**3, RKIND)
                 CO2(i,j,k)     = CO2(i,j,k)*real(aye**3, RKIND)
                 OI(i,j,k)      = OI(i,j,k)*real(aye**3, RKIND)
                 OH(i,j,k)      = OH(i,j,k)*real(aye**3, RKIND)
                 H2O(i,j,k)     = H2O(i,j,k)*real(aye**3, RKIND)
                 O2(i,j,k)      = O2(i,j,k)*real(aye**3, RKIND)
                 SiI(i,j,k)     = SiI(i,j,k)*real(aye**3, RKIND)
                 SiOI(i,j,k)    = SiOI(i,j,k)*real(aye**3, RKIND)
                 SiO2I(i,j,k)   = SiO2I(i,j,k)*real(aye**3, RKIND)
                 CH(i,j,k)      = CH(i,j,k)*real(aye**3, RKIND)
                 CH2(i,j,k)     = CH2(i,j,k)*real(aye**3, RKIND)
                 COII(i,j,k)    = COII(i,j,k)*real(aye**3, RKIND)
                 OII(i,j,k)     = OII(i,j,k)*real(aye**3, RKIND)
                 OHII(i,j,k)    = OHII(i,j,k)*real(aye**3, RKIND)
                 H2OII(i,j,k)   = H2OII(i,j,k)*real(aye**3, RKIND)
                 H3OII(i,j,k)   = H3OII(i,j,k)*real(aye**3, RKIND)
                 O2II(i,j,k)    = O2II(i,j,k)*real(aye**3, RKIND)
                 if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
                 if (idspecies .gt. 0) then
                 Mg(i,j,k)      = Mg(i,j,k)*real(aye**3, RKIND)
                 endif
                 if (idspecies .gt. 1) then
                 Al(i,j,k)      = Al(i,j,k)*real(aye**3, RKIND)
                 S(i,j,k)       = S(i,j,k)*real(aye**3, RKIND)
                 Fe(i,j,k)      = Fe(i,j,k)*real(aye**3, RKIND)
                 endif
                 endif
              endif
!!            endif
           enddo
        endif
        if (idustfield .eq. 1) then
           do i = is+1, ie+1
              dust(i,j,k) = dust(i,j,k)*real(aye**3, RKIND)
              if ( ( igrgr .eq. 1 ) .or. ( idsub .eq. 1 ) ) then
!!            if (metal(i,j,k) .gt. 1.d-9 * d(i,j,k)) then
              if (idspecies .gt. 0) then
                 MgSiO3(i,j,k)  = MgSiO3(i,j,k)*real(aye**3, RKIND)
                 AC(i,j,k)      = AC(i,j,k)*real(aye**3, RKIND)
              endif
              if (idspecies .gt. 1) then
                 SiM(i,j,k)     = SiM(i,j,k)*real(aye**3, RKIND)
                 FeM(i,j,k)     = FeM(i,j,k)*real(aye**3, RKIND)
                 Mg2SiO4(i,j,k) = Mg2SiO4(i,j,k)*real(aye**3, RKIND)
                 Fe3O4(i,j,k)   = Fe3O4(i,j,k)*real(aye**3, RKIND)
                 SiO2D(i,j,k)   = SiO2D(i,j,k)*real(aye**3, RKIND)
                 MgO(i,j,k)     = MgO(i,j,k)*real(aye**3, RKIND)
                 FeS(i,j,k)     = FeS(i,j,k)*real(aye**3, RKIND)
                 Al2O3(i,j,k)   = Al2O3(i,j,k)*real(aye**3, RKIND)
              endif
              if (idspecies .gt. 2) then
                 reforg(i,j,k)  = reforg(i,j,k)*real(aye**3, RKIND)
                 volorg(i,j,k)  = volorg(i,j,k)*real(aye**3, RKIND)
                 H2Oice(i,j,k)  = H2Oice(i,j,k)*real(aye**3, RKIND)
              endif
!!            endif
              endif
           enddo
        endif
      enddo
#ifdef _OPENMP
!$omp end parallel do
#endif

      endif

      return
      end
