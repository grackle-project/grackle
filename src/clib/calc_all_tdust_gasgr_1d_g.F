      subroutine calc_all_tdust_gasgr_1d_g(in, jn, kn, nratec,
     &                idustfield, is, ie, j, k, fgr, gamma_isrfa,
     &                trad, gasgra, indixe, tdef, tgas, tdust,
     &                metallicity, dust2gas, nh, gasgr_tdust,
     &                itmask_metal,
     &                idspecies, itdmulti, gr_N, gr_Size, gr_dT,
     &                gr_Td, tSiM, tFeM, tMg2SiO4, tMgSiO3, tFe3O4,
     &                tAC, tSiO2D, tMgO, tFeS, tAl2O3, treforg,
     &                tvolorg, tH2Oice, gasgr2a, gamma_isrf2a,
     &                coolunit, gasgr, myisrf, sgSiM, sgFeM, sgMg2SiO4,
     &                sgMgSiO3, sgFe3O4, sgAC, sgSiO2D, sgMgO, sgFeS,
     &                sgAl2O3, sgreforg, sgvolorg, sgH2Oice, sgtot,
     &                alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC,
     &                alSiO2D, alMgO, alFeS, alAl2O3, alreforg,
     &                alvolorg, alH2Oice, altot, kpSiM, kpFeM,
     &                kpMg2SiO4, kpMgSiO3, kpFe3O4, kpAC, kpSiO2D,
     &                kpMgO, kpFeS, kpAl2O3, kpreforg, kpvolorg,
     &                kpH2Oice, kptot, gasSiM, gasFeM, gasMg2SiO4,
     &                gasMgSiO3, gasFe3O4, gasAC, gasSiO2D, gasMgO,
     &                gasFeS, gasAl2O3, gasreforg, gasvolorg,
     &                gasH2Oice
     &              )

!  PURPOSE:
!    Calculate all dust temperature(s) and the gas to grain heat
!    transfer rate(s) (the latter is commonly called gasgr)
!
!    An argument could be made for directly using gasgr to compute
!    contributions to edot within this routine, rather than returning
!    the gasgr value(s). But that should be reconsidered in the future.
!
!    We could significantly reduce the amount of buffer space allocated
!    within the routine. But, we will hold off on doing that until after
!    we have transcribed to C/C++
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!

!-----------------------------------------------------------------------

      implicit NONE
#include "grackle_fortran_types.def"
#include "phys_const.def"

!  Arguments

      integer in, jn, kn, is, ie, j, k, nratec, idustfield

      real*8 fgr
      real*8 gasgra(nratec), gamma_isrfa

!  Parameters

      real*8, parameter :: mh = mass_h      !DPC

!  Locals

      integer i
      real*8 trad, coolunit

!  Slice locals
 
      integer*8 indixe(in)
      real*8 tdef(in), tgas(in), tdust(in), nh(in), metallicity(in),
     &       dust2gas(in)

!  Cooling/heating slice locals

      real*8 gasgr(in), gasgr_tdust(in), myisrf(in)
      integer idspecies, itdmulti
!     opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
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
!  Iteration mask

      MASK_TYPE itmask_metal(in)

          ! Calculate heating from interstellar radiation field
          !  -> this is ONLY used when `itmask_metal .eq. MASK_TRUE`

      do i = is+1, ie+1
         if ( itmask_metal(i) .ne. MASK_FALSE ) then

         if (idspecies .eq. 0 ) then
            if (idustfield .gt. 0) then
               mygisrf(i) = gamma_isrfa
     &                    * fgr / dust2gas(i) * metallicity(i)
               !! correct with the depletion or enhancement of condensation rate.
            else
               mygisrf(i) = gamma_isrfa
            endif

         else ! idspecies

            if (itdmulti .eq. 0) then
            
               mygisrf(i) = gamma_isrf2a * sgtot(i)
               !! in UV, absorption coefficient Q ~ 1 (Goldsmith 2001)
               !! so we use the geometrical cross-section of grains [cgs]

            else

            if (idspecies .gt. 0) then
               gisrfMgSiO3  (i) = gamma_isrf2a * sgMgSiO3  (i)
!!             write(*,*) 'sil', d(i,j,k), gamma_isrf2a, sgMgSiO3(i)
               gisrfAC      (i) = gamma_isrf2a * sgAC      (i)
!!             write(*,*) 'car', d(i,j,k), gamma_isrf2a, sgMgSiO3(i)
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

         endif ! idspecies

         endif
      enddo

!     --- Gas to grain heat transfer ---

!     Look up gas/grain heat transfer rates

         do i = is+1, ie+1
            if ( itmask_metal(i) .ne. MASK_FALSE ) then

            if(idspecies .eq. 0) then

               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))

!!             gasgr_tdust(i) = fgr * gasgr(i) * coolunit / mh
               gasgr_tdust(i) = (dust2gas(i) / metallicity(i))
     &                              * gasgr(i) * coolunit / mh
                    !! apply to (idustfield .eq. 1) GC20200701

            else ! idspecies

               fv2k    = gasgr2a(indixe(i)) + tdef(i)
     &              *(gasgr2a(indixe(i)+1) -gasgr2a(indixe(i)))

               fac = coolunit / mh 

               if ( itdmulti .eq. 0 ) then

               gasgr(i) = fv2k * sgtot(i)

               gasgr_tdust(i) = gasgr(i) * fac
                 
               else

               if (idspecies .gt. 0) then
               gasMgSiO3  (i) = fv2k * sgMgSiO3  (i)
               gasAC      (i) = fv2k * sgAC      (i)
               endif
               if (idspecies .gt. 1) then
               gasSiM     (i) = fv2k * sgSiM     (i)
               gasFeM     (i) = fv2k * sgFeM     (i)
               gasMg2SiO4 (i) = fv2k * sgMg2SiO4 (i)
               gasFe3O4   (i) = fv2k * sgFe3O4   (i)
               gasSiO2D   (i) = fv2k * sgSiO2D   (i)
               gasMgO     (i) = fv2k * sgMgO     (i)
               gasFeS     (i) = fv2k * sgFeS     (i)
               gasAl2O3   (i) = fv2k * sgAl2O3   (i)
               endif
               if (idspecies .gt. 2) then
               gasreforg  (i) = fv2k * sgreforg  (i)
               gasvolorg  (i) = fv2k * sgvolorg  (i)
               gasH2Oice  (i) = fv2k * sgH2Oice  (i)
               endif

               if (idspecies .gt. 0) then
               gasgr_tMgSiO3  (i) = gasMgSiO3  (i) * fac
               gasgr_tAC      (i) = gasAC      (i) * fac
               endif
               if (idspecies .gt. 1) then
               gasgr_tSiM     (i) = gasSiM     (i) * fac
               gasgr_tFeM     (i) = gasFeM     (i) * fac
               gasgr_tMg2SiO4 (i) = gasMg2SiO4 (i) * fac
               gasgr_tFe3O4   (i) = gasFe3O4   (i) * fac
               gasgr_tSiO2D   (i) = gasSiO2D   (i) * fac
               gasgr_tMgO     (i) = gasMgO     (i) * fac
               gasgr_tFeS     (i) = gasFeS     (i) * fac
               gasgr_tAl2O3   (i) = gasAl2O3   (i) * fac
               endif
               if (idspecies .gt. 2) then
               gasgr_treforg  (i) = gasreforg  (i) * fac
               gasgr_tvolorg  (i) = gasvolorg  (i) * fac
               gasgr_tH2Oice  (i) = gasH2Oice  (i) * fac
               endif

               endif

            endif ! idspecies

            endif
         enddo

!     Compute dust temperature

         if (itdmulti .eq. 0) then

            call calc_tdust_1d_g(tdust, tgas, nh, gasgr_tdust,
     &           mygisrf, myisrf, itmask_metal, trad, in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, altot, kptot, idspecies)

         else

            if (idspecies .gt. 0) then
            call calc_tdust_1d_g(tMgSiO3  , tgas, nh, gasgr_tMgSiO3  ,
     &           gisrfMgSiO3  , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgSiO3  , kpMgSiO3   
     &         , idspecies)
          
            call calc_tdust_1d_g(tAC      , tgas, nh, gasgr_tAC      ,
     &           gisrfAC      , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAC      , kpAC       
     &         , idspecies)
            endif

            if (idspecies .gt. 1) then
            call calc_tdust_1d_g(tSiM     , tgas, nh, gasgr_tSiM     ,
     &           gisrfSiM     , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiM     , kpSiM      
     &         , idspecies)
          
            call calc_tdust_1d_g(tFeM     , tgas, nh, gasgr_tFeM     ,
     &           gisrfFeM     , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeM     , kpFeM      
     &         , idspecies)
          
            call calc_tdust_1d_g(tMg2SiO4 , tgas, nh, gasgr_tMg2SiO4 ,
     &           gisrfMg2SiO4 , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMg2SiO4 , kpMg2SiO4  
     &         , idspecies)
          
            call calc_tdust_1d_g(tFe3O4   , tgas, nh, gasgr_tFe3O4   ,
     &           gisrfFe3O4   , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFe3O4   , kpFe3O4    
     &         , idspecies)
          
            call calc_tdust_1d_g(tSiO2D   , tgas, nh, gasgr_tSiO2D   ,
     &           gisrfSiO2D   , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alSiO2D   , kpSiO2D    
     &         , idspecies)
          
            call calc_tdust_1d_g(tMgO     , tgas, nh, gasgr_tMgO     ,
     &           gisrfMgO     , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alMgO     , kpMgO      
     &         , idspecies)
          
            call calc_tdust_1d_g(tFeS     , tgas, nh, gasgr_tFeS     ,
     &           gisrfFeS     , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alFeS     , kpFeS      
     &         , idspecies)
          
            call calc_tdust_1d_g(tAl2O3   , tgas, nh, gasgr_tAl2O3   ,
     &           gisrfAl2O3   , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alAl2O3   , kpAl2O3    
     &         , idspecies)
            endif

            if (idspecies .gt. 2) then
            call calc_tdust_1d_g(treforg  , tgas, nh, gasgr_treforg  ,
     &           gisrfreforg  , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alreforg  , kpreforg   
     &         , idspecies)
          
            call calc_tdust_1d_g(tvolorg  , tgas, nh, gasgr_tvolorg  ,
     &           gisrfvolorg  , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alvolorg  , kpvolorg   
     &         , idspecies)
          
            call calc_tdust_1d_g(tH2Oice  , tgas, nh, gasgr_tH2Oice  ,
     &           gisrfH2Oice  , myisrf, itmask_metal, trad, 
     &           in, is, ie, j, k
     &         , gr_N, gr_Size, gr_dT, gr_Td, alH2Oice  , kpH2Oice   
     &         , idspecies)
            endif

         endif
      end
