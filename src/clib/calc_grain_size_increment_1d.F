#include "phys_const.def"

!     Compute grain size increment

      subroutine calc_grain_size_increment_1d(
     &                immulti, imabund, idspecies, igrgr, itmask
     &              , in, jn, kn, is, ie, j, k, dom, d
     &              , SiM, FeM, Mg2SiO4, MgSiO3, Fe3O4
     &              , AC, SiO2D, MgO, FeS, Al2O3
     &              , reforg, volorg, H2Oice
     &              , metal, metal_loc
     &              , metal_C13, metal_C20, metal_C25, metal_C30
     &              , metal_F13, metal_F15, metal_F50, metal_F80
     &              , metal_P170, metal_P200, metal_Y19
     &              , SN0_N
     &              , SN0_fSiM, SN0_fFeM, SN0_fMg2SiO4, SN0_fMgSiO3
     &              , SN0_fFe3O4, SN0_fAC, SN0_fSiO2D, SN0_fMgO
     &              , SN0_fFeS, SN0_fAl2O3
     &              , SN0_freforg, SN0_fvolorg, SN0_fH2Oice
     &              , SN0_r0SiM, SN0_r0FeM, SN0_r0Mg2SiO4, SN0_r0MgSiO3
     &              , SN0_r0Fe3O4, SN0_r0AC, SN0_r0SiO2D, SN0_r0MgO
     &              , SN0_r0FeS, SN0_r0Al2O3
     &              , SN0_r0reforg, SN0_r0volorg, SN0_r0H2Oice
     &              , gr_N, gr_Size, gr_dT, gr_Td
     &              , SN0_kpSiM, SN0_kpFeM, SN0_kpMg2SiO4, SN0_kpMgSiO3
     &              , SN0_kpFe3O4, SN0_kpAC, SN0_kpSiO2D, SN0_kpMgO
     &              , SN0_kpFeS, SN0_kpAl2O3
     &              , SN0_kpreforg, SN0_kpvolorg, SN0_kpH2Oice
     &              , sgSiM, sgFeM, sgMg2SiO4, sgMgSiO3, sgFe3O4, sgAC
     &              , sgSiO2D, sgMgO, sgFeS, sgAl2O3
     &              , sgreforg, sgvolorg, sgH2Oice, sgtot
     &              , alSiM, alFeM, alMg2SiO4, alMgSiO3, alFe3O4, alAC
     &              , alSiO2D, alMgO, alFeS, alAl2O3
     &              , alreforg, alvolorg, alH2Oice, altot
     &            )

      implicit NONE
#include "grackle_fortran_types.def"
#include "dust_const.def"

!     in
      integer in, jn, kn, is, ie, j, k
      MASK_TYPE itmask(in)
      integer immulti, imabund, idspecies, igrgr
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  SiM(in,jn,kn), FeM(in,jn,kn), Mg2SiO4(in,jn,kn)
     &      , MgSiO3(in,jn,kn), Fe3O4(in,jn,kn), AC(in,jn,kn)
     &      , SiO2D(in,jn,kn), MgO(in,jn,kn), FeS(in,jn,kn)
     &      , Al2O3(in,jn,kn)
     &      , reforg(in,jn,kn), volorg(in,jn,kn), H2Oice(in,jn,kn)
      R_PREC  metal(in,jn,kn)
     &      , metal_loc(in,jn,kn)
     &      , metal_C13(in,jn,kn), metal_C20(in,jn,kn)
     &      , metal_C25(in,jn,kn), metal_C30(in,jn,kn)
     &      , metal_F13(in,jn,kn), metal_F15(in,jn,kn)
     &      , metal_F50(in,jn,kn), metal_F80(in,jn,kn)
     &      , metal_P170(in,jn,kn), metal_P200(in,jn,kn)
     &      , metal_Y19(in,jn,kn)
!     table
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
!     out
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
!     array
      integer iSN, nSN, iSN0
      integer SN_i(SN0_N)
      R_PREC  SN_metal(in, SN0_N)
      real*8  SN_fSiM(SN0_N), SN_fFeM(SN0_N)
     &      , SN_fMg2SiO4(SN0_N), SN_fMgSiO3(SN0_N)
     &      , SN_fFe3O4(SN0_N), SN_fAC(SN0_N)
     &      , SN_fSiO2D(SN0_N), SN_fMgO(SN0_N)
     &      , SN_fFeS(SN0_N), SN_fAl2O3(SN0_N)
     &      , SN_freforg(SN0_N)
     &      , SN_fvolorg(SN0_N), SN_fH2Oice(SN0_N)
      real*8  SN_r0SiM(3,SN0_N), SN_r0FeM(3,SN0_N)
     &      , SN_r0Mg2SiO4(3,SN0_N), SN_r0MgSiO3(3,SN0_N)
     &      , SN_r0Fe3O4(3,SN0_N), SN_r0AC(3,SN0_N)
     &      , SN_r0SiO2D(3,SN0_N), SN_r0MgO(3,SN0_N)
     &      , SN_r0FeS(3,SN0_N), SN_r0Al2O3(3,SN0_N)
     &      , SN_r0reforg(3,SN0_N)
     &      , SN_r0volorg(3,SN0_N), SN_r0H2Oice(3,SN0_N)
      real*8  SN_kpSiM(gr_Size,SN0_N), SN_kpFeM(gr_Size,SN0_N)
     &      , SN_kpMg2SiO4(gr_Size,SN0_N), SN_kpMgSiO3(gr_Size,SN0_N)
     &      , SN_kpFe3O4(gr_Size,SN0_N), SN_kpAC(gr_Size,SN0_N)
     &      , SN_kpSiO2D(gr_Size,SN0_N), SN_kpMgO(gr_Size,SN0_N)
     &      , SN_kpFeS(gr_Size,SN0_N), SN_kpAl2O3(gr_Size,SN0_N)
     &      , SN_kpreforg(gr_Size,SN0_N)
     &      , SN_kpvolorg(gr_Size,SN0_N), SN_kpH2Oice(gr_Size,SN0_N)
!     local
      integer i
!!    integer iTd, iTd0

!!    iSN0 = 2
!!    write(*,*) SN0_fMgSiO3(iSN0)
!!    write(*,*) SN0_r0MgSiO3  (:,iSN0)
!!    do iTd = 1, gr_N(2)
!!       iTd0 = (iTd - 1) * gr_N(1)
!!       write(*,*) SN0_kpMgSiO3  (iTd0+1,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+2,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+3,iSN0)
!!   &            , SN0_kpMgSiO3  (iTd0+4,iSN0)
!!    enddo
!!    stop
        
      ! make arrays
      if (immulti .eq. 0) then
         nSN = 1
         SN_i(1) = imabund + 1
         SN_metal(:,1) = metal(:,j,k)
      else
         nSN = 0
         if(maxval(metal_loc(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  1
            SN_metal(:,nSN) = metal_loc(:,j,k)
         endif
         if(maxval(metal_C13(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  2
            SN_metal(:,nSN) = metal_C13(:,j,k)
         endif
         if(maxval(metal_C20(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  3
            SN_metal(:,nSN) = metal_C20(:,j,k)
         endif
         if(maxval(metal_C25(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  4
            SN_metal(:,nSN) = metal_C25(:,j,k)
         endif
         if(maxval(metal_C30(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  5
            SN_metal(:,nSN) = metal_C30(:,j,k)
         endif
         if(maxval(metal_F13(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  6
            SN_metal(:,nSN) = metal_F13(:,j,k)
         endif
         if(maxval(metal_F15(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  7
            SN_metal(:,nSN) = metal_F15(:,j,k)
         endif
         if(maxval(metal_F50(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  8
            SN_metal(:,nSN) = metal_F50(:,j,k)
         endif
         if(maxval(metal_F80(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) =  9
            SN_metal(:,nSN) = metal_F80(:,j,k)
         endif
         if(maxval(metal_P170(:,j,k)/metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 10
            SN_metal(:,nSN) = metal_P170(:,j,k)
         endif
         if(maxval(metal_P200(:,j,k)/metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 11
            SN_metal(:,nSN) = metal_P200(:,j,k)
         endif
         if(maxval(metal_Y19(:,j,k) /metal(:,j,k)) .gt. 0.01_DKIND) then
            nSN = nSN + 1
            SN_i(nSN) = 12
            SN_metal(:,nSN) = metal_Y19(:,j,k)
         endif
      endif
         
      do iSN = 1, nSN
         iSN0 = SN_i(iSN)
         if ( idspecies .gt. 0 ) then
            SN_fMgSiO3     (iSN) = SN0_fMgSiO3     (iSN0)
            SN_fAC         (iSN) = SN0_fAC         (iSN0)
            SN_r0MgSiO3  (:,iSN) = SN0_r0MgSiO3  (:,iSN0)
            SN_r0AC      (:,iSN) = SN0_r0AC      (:,iSN0)
            SN_kpMgSiO3  (:,iSN) = SN0_kpMgSiO3  (:,iSN0)
            SN_kpAC      (:,iSN) = SN0_kpAC      (:,iSN0)
         endif
         if ( idspecies .gt. 1 ) then
            SN_fSiM        (iSN) = SN0_fSiM        (iSN0)
            SN_fFeM        (iSN) = SN0_fFeM        (iSN0)
            SN_fMg2SiO4    (iSN) = SN0_fMg2SiO4    (iSN0)
            SN_fFe3O4      (iSN) = SN0_fFe3O4      (iSN0)
            SN_fSiO2D      (iSN) = SN0_fSiO2D      (iSN0)
            SN_fMgO        (iSN) = SN0_fMgO        (iSN0)
            SN_fFeS        (iSN) = SN0_fFeS        (iSN0)
            SN_fAl2O3      (iSN) = SN0_fAl2O3      (iSN0)
            SN_r0SiM     (:,iSN) = SN0_r0SiM     (:,iSN0)
            SN_r0FeM     (:,iSN) = SN0_r0FeM     (:,iSN0)
            SN_r0Mg2SiO4 (:,iSN) = SN0_r0Mg2SiO4 (:,iSN0)
            SN_r0Fe3O4   (:,iSN) = SN0_r0Fe3O4   (:,iSN0)
            SN_r0SiO2D   (:,iSN) = SN0_r0SiO2D   (:,iSN0)
            SN_r0MgO     (:,iSN) = SN0_r0MgO     (:,iSN0)
            SN_r0FeS     (:,iSN) = SN0_r0FeS     (:,iSN0)
            SN_r0Al2O3   (:,iSN) = SN0_r0Al2O3   (:,iSN0)
            SN_kpSiM     (:,iSN) = SN0_kpSiM     (:,iSN0)
            SN_kpFeM     (:,iSN) = SN0_kpFeM     (:,iSN0)
            SN_kpMg2SiO4 (:,iSN) = SN0_kpMg2SiO4 (:,iSN0)
            SN_kpFe3O4   (:,iSN) = SN0_kpFe3O4   (:,iSN0)
            SN_kpSiO2D   (:,iSN) = SN0_kpSiO2D   (:,iSN0)
            SN_kpMgO     (:,iSN) = SN0_kpMgO     (:,iSN0)
            SN_kpFeS     (:,iSN) = SN0_kpFeS     (:,iSN0)
            SN_kpAl2O3   (:,iSN) = SN0_kpAl2O3   (:,iSN0)
         endif
         if ( idspecies .gt. 2 ) then
            SN_freforg     (iSN) = SN0_freforg     (iSN0)
            SN_fvolorg     (iSN) = SN0_fvolorg     (iSN0)
            SN_fH2Oice     (iSN) = SN0_fH2Oice     (iSN0)
            SN_r0reforg  (:,iSN) = SN0_r0reforg  (:,iSN0)
            SN_r0volorg  (:,iSN) = SN0_r0volorg  (:,iSN0)
            SN_r0H2Oice  (:,iSN) = SN0_r0H2Oice  (:,iSN0)
            SN_kpreforg  (:,iSN) = SN0_kpreforg  (:,iSN0)
            SN_kpvolorg  (:,iSN) = SN0_kpvolorg  (:,iSN0)
            SN_kpH2Oice  (:,iSN) = SN0_kpH2Oice  (:,iSN0)
         endif
      enddo

      !! calculate size increment

      if (idspecies .gt. 0) then
!!    write(*,*) 'MgSiO3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgSiO3  , SN_metal, SN_fMgSiO3  , SN_r0MgSiO3  
     &        , sMgSiO3  , sgMgSiO3  , alMgSiO3  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgSiO3
     &            )

!!    write(*,*) 'AC'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, AC      , SN_metal, SN_fAC      , SN_r0AC      
     &        , sAC      , sgAC      , alAC      
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAC
     &            )
      endif

      if (idspecies .gt. 1) then
!!    write(*,*) 'SiM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiM     , SN_metal, SN_fSiM     , SN_r0SiM     
     &        , sSiM     , sgSiM     , alSiM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiM
     &            )

!!    write(*,*) 'FeM'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeM     , SN_metal, SN_fFeM     , SN_r0FeM     
     &        , sFeM     , sgFeM     , alFeM     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeM
     &            )

!!    write(*,*) 'Mg2SiO4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Mg2SiO4 , SN_metal, SN_fMg2SiO4 , SN_r0Mg2SiO4 
     &        , sMg2SiO4 , sgMg2SiO4 , alMg2SiO4 
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMg2SiO4
     &            )

!!    write(*,*) 'Fe3O4'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Fe3O4   , SN_metal, SN_fFe3O4   , SN_r0Fe3O4   
     &        , sFe3O4   , sgFe3O4   , alFe3O4   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFe3O4
     &            )

!!    write(*,*) 'SiO2D'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, SiO2D   , SN_metal, SN_fSiO2D   , SN_r0SiO2D   
     &        , sSiO2D   , sgSiO2D   , alSiO2D   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpSiO2D
     &            )

!!    write(*,*) 'MgO'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, MgO     , SN_metal, SN_fMgO     , SN_r0MgO     
     &        , sMgO     , sgMgO     , alMgO     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpMgO
     &            )

!!    write(*,*) 'FeS'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, FeS     , SN_metal, SN_fFeS     , SN_r0FeS     
     &        , sFeS     , sgFeS     , alFeS     
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpFeS
     &            )

!!    write(*,*) 'Al2O3'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, Al2O3   , SN_metal, SN_fAl2O3   , SN_r0Al2O3   
     &        , sAl2O3   , sgAl2O3   , alAl2O3   
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpAl2O3
     &            )
      endif

      if (idspecies .gt. 2) then
!!    write(*,*) 'reforg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, reforg  , SN_metal, SN_freforg  , SN_r0reforg  
     &        , sreforg  , sgreforg  , alreforg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpreforg
     &            )

!!    write(*,*) 'volorg'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, volorg  , SN_metal, SN_fvolorg  , SN_r0volorg  
     &        , svolorg  , sgvolorg  , alvolorg  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpvolorg
     &            )

!!    write(*,*) 'H2Oice'
      call calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, H2Oice  , SN_metal, SN_fH2Oice  , SN_r0H2Oice  
     &        , sH2Oice  , sgH2Oice  , alH2Oice  
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kpH2Oice
     &            )
      endif

      do i = is+1, ie+1
      if ( itmask(i) .ne. MASK_FALSE ) then
            
         if (idspecies .gt. 0) then
            sgtot  (i) = sgMgSiO3    (i)
     &                 + sgAC        (i)
         endif
         if (idspecies .gt. 1) then
            sgtot  (i) = sgtot       (i)
     &                 + sgSiM       (i)
     &                 + sgFeM       (i)
     &                 + sgMg2SiO4   (i)
     &                 + sgFe3O4     (i)
     &                 + sgSiO2D     (i)
     &                 + sgMgO       (i)
     &                 + sgFeS       (i)
     &                 + sgAl2O3     (i)
         endif
         if (idspecies .gt. 2) then
            sgtot  (i) = sgtot(i)
     &                 + sgreforg    (i)
     &                 + sgvolorg    (i)
     &                 + sgH2Oice    (i)
         endif
            
         if (idspecies .gt. 0) then
            altot(:,i) = alMgSiO3  (:,i)
     &                 + alAC      (:,i)
         endif
         if (idspecies .gt. 1) then
            altot(:,i) = altot     (:,i)
     &                 + alSiM     (:,i)
     &                 + alFeM     (:,i)
     &                 + alMg2SiO4 (:,i)
     &                 + alFe3O4   (:,i)
     &                 + alSiO2D   (:,i)
     &                 + alMgO     (:,i)
     &                 + alFeS     (:,i)
     &                 + alAl2O3   (:,i)
         endif
         if (idspecies .gt. 2) then
            altot(:,i) = altot     (:,i)
     &                 + alreforg  (:,i)
     &                 + alvolorg  (:,i)
     &                 + alH2Oice  (:,i)
         endif

      endif
      enddo

      return
      end

!=======================================================================
!//////////////  CALC_GRAIN_SIZE_INCREMENT_SPECIES_1D \\\\\\\\\\\\\\\\\\

      subroutine calc_grain_size_increment_species_1d(
     &          igrgr, itmask, SN0_N
     &        , in, jn, kn, is, ie, j, k, dom, d
     &        , nSN, dsp, SN_metal, SN_fsp, SN_r0sp
     &        , ssp, sgsp, alsp
     &        , gr_N, gr_Size, gr_dT, gr_Td, SN_kp0sp
     &            )


      implicit NONE

#include "grackle_fortran_types.def"

! input
      integer in, jn, kn, is, ie, j, k
      MASK_TYPE itmask(in)
      integer igrgr
      integer iSN, nSN, SN0_N
      real*8  dom
      R_PREC  d(in,jn,kn)
      R_PREC  dsp(in,jn,kn), SN_metal(in, SN0_N)
! table
      real*8  SN_fsp(SN0_N), SN_r0sp(3,SN0_N)
      real*8  ssp
! opacity table
      integer gr_N(2), gr_Size
      real*8  gr_dT, gr_Td(gr_N(2))
      real*8  SN_kp0sp(gr_Size,SN0_N)
! output
      real*8 sgsp(in), alsp(gr_N(2),in)
! local
      integer i
      real*8 coef0, coef1, coef2, coef3
      real*8 dsp0, SN_sgsp, SN_kpsp
      real*8 SN_dsp0(SN0_N), SN_nsp0(SN0_N)
      real*8 drsp(in)
      real*8, parameter :: pi = pi_val
      real*8, parameter :: mh = mass_h
!     debug
      real*8 SN_dsp(SN0_N), SN_msp(SN0_N), dsp1
      integer iTd, iTd0

!       write(*,*) SN_fsp(1)
!       write(*,*) SN_r0sp(1,1)
!    &           , SN_r0sp(2,1)
!    &           , SN_r0sp(3,1)
!     do iTd = 1, gr_N(2)
!       iTd0 = (iTd-1)*gr_N(1)
!       write(*,*) gr_Td(iTd)
!    &           , SN_kp0sp(iTd0+1,1)
!    &           , SN_kp0sp(iTd0+2,1)
!    &           , SN_kp0sp(iTd0+3,1)
!    &           , SN_kp0sp(iTd0+4,1)
!     enddo
 
      do i = is+1, ie+1
         if ( itmask(i) .ne. MASK_FALSE ) then
!!!!!!!!!!!!!
!!!!     if( dsp(i,j,k) .gt. 1.d-15*d(i,j,k) ) then
!!!!!!!!!!!!!
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_dsp0(iSN) = SN_fsp(iSN) * SN_metal(i, iSN)
!!                write(*,*)  iSN, SN_fsp(iSN)
!!   &                      , SN_metal(i, iSN) 
               endif
            enddo

            if(igrgr .eq. 0) then

            drsp(i) = 0.d0

            else

!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!!!! TEST !!!
!!!         dsp(i,j,k) = 3.d0 * dsp(i,j,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculate grain size increment from the conservation of mass

            coef0 = 0.d0
            coef1 = 0.d0
            coef2 = 0.d0
            coef3 = 0.d0
            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN) / SN_r0sp(3,iSN)
                  coef0 = coef0 +        SN_dsp0(iSN)
                  coef1 = coef1 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN)
                  coef2 = coef2 + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN)
                  coef3 = coef3 +        SN_nsp0(iSN)
               endif
            enddo
!!          write(*,*) coef0, dsp(i,j,k)
            coef0 = coef0 - dsp(i,j,k)
!!          write(*,*) coef0, coef1, coef2, coef3

            coef0 = coef0 / coef3
            coef1 = coef1 / coef3
            coef2 = coef2 / coef3
            
            call solve_cubic_equation(coef2, coef1, coef0, drsp(i))
!!          write(*,*) coef0, coef1, coef2, coef3, drsp(i)

            drsp(i) = max(drsp(i), 0.d0)

            endif


! calculate number density (code_density / g)

            do iSN = 1, nSN
               if(SN_fsp(iSN) .gt. 0.d0) then
                  SN_nsp0(iSN) = SN_dsp0(iSN)
     &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
!                 if(iSN.eq.1)
!    &              write(*,*) 'aa', SN_dsp0(iSN), dsp(i,j,k)
!    &             , (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
!    &                      , SN_nsp0(iSN)
               else
                  SN_nsp0(iSN) = 0.d0
               endif
            enddo

! calculate geometrical cross-section per unit gas mass
            sgsp(i) = 0.d0
            do iSN = 1, nSN
               if( SN_fsp(iSN) .gt. 0.d0) then
                  SN_sgsp = pi *
     &               (        SN_r0sp(2,iSN)
     &               + 2.d0 * SN_r0sp(1,iSN) * drsp(i)
     &               +                         drsp(i)**2
     &               )
               else
                  SN_sgsp = 0.d0
               endif
               sgsp(i) = sgsp(i) + SN_nsp0(iSN) * SN_sgsp
!!             if(iSN.eq.1)
!!   &           write(*,*) dsp(i,j,k)/d(i,j,k)
!!   &                  , SN_dsp0(iSN)/d(i,j,k)
!!   &              , ssp, SN_r0sp(3, iSN), SN_r0sp(2,iSN) 
!!   &                  , drsp(i), sgsp(i)/d(i,j,k)
            enddo
            sgsp(i) = sgsp(i) / d(i,j,k)
    
! calculate optical cross-section per unit gas mass
            do iTd = 1, gr_N(2)
               iTd0 = (iTd-1)*gr_N(1)
               alsp(iTd,i) = 0.d0
               do iSN = 1, nSN
                  if( SN_fsp(iSN) .gt. 0.d0) then
                     SN_kpsp = 4.d0*pi/3.d0 * ssp *
     &                  (        SN_kp0sp(iTd0+4,iSN)
     &                  + 3.d0 * SN_kp0sp(iTd0+3,iSN) * drsp(i)
     &                  + 3.d0 * SN_kp0sp(iTd0+2,iSN) * drsp(i)**2
     &                  +        SN_kp0sp(iTd0+1,iSN) * drsp(i)**3
     &                  )
                  else
                     SN_kpsp = 0.d0
                  endif
                  alsp(iTd,i) = alsp(iTd,i) + SN_nsp0(iSN) * SN_kpsp
!                 if(iSN.eq.1)
!    &               write(*,*) gr_Td(iTd)
!    &                        , SN_dsp0(iSN)/d(i,j,k)
!    &             , (4.d0*pi/3.d0 * ssp * SN_kp0sp(iTd0+4,iSN))
!    &             / (4.d0*pi/3.d0 * ssp * SN_r0sp(3,iSN))
               enddo
               alsp(iTd,i) = alsp(iTd,i) / d(i,j,k)
            enddo

! decompose grain mass density for each SN model !!!!!!!!!!!!!!!!!!!!!!!
!           dsp0 = 0.d0
!           dsp1 = 0.d0
!           do iSN = 1, nSN
!              dsp0 = dsp0 + SN_dsp0(iSN)
!              if( SN_fsp(iSN) .gt. 0.d0) then
!                 SN_dsp(iSN) = 4.d0*pi/3.d0 * ssp * (
!    &                        SN_nsp0(iSN) * SN_r0sp(3,iSN)
!    &               + 3.d0 * SN_nsp0(iSN) * SN_r0sp(2,iSN) * drsp(i)
!    &               + 3.d0 * SN_nsp0(iSN) * SN_r0sp(1,iSN) * drsp(i)**2
!    &               +        SN_nsp0(iSN)                  * drsp(i)**3
!    &               )
!                 SN_msp(iSN) = 4.d0*pi/3.d0 * ssp * (
!    &                        SN_r0sp(3,iSN)
!    &               + 3.d0 * SN_r0sp(2,iSN) * drsp(i)
!    &               + 3.d0 * SN_r0sp(1,iSN) * drsp(i)**2
!    &               +                         drsp(i)**3
!    &               )
!!                if(iSN.eq.2)
!!   &            write(*,*) iSN, (SN_dsp(iSN)/SN_msp(iSN))/SN_nsp0(iSN)
!                 dsp1 = dsp1 + SN_dsp(iSN)
!              endif
!           enddo
!           write(*,*) drsp(i), dsp1/dsp(i,j,k)
!!!!!!!!!!!!!
!!!!     else

!!!!        sgsp  (i) = tiny8
!!!!        alsp(:,i) = tiny8

!!!!     endif
!!!!!!!!!!!!!

         endif ! itmask
      enddo

      return
      end

!=======================================================================
!////////////////  SUBROUTINE SOLVE_CUBIC_EQUATION \\\\\\\\\\\\\\\\\\\\\

      subroutine solve_cubic_equation(a, b, c, root)

      implicit NONE
#include "grackle_fortran_types.def"

!  General Arguments

      real*8 a, b, c
      real*8 root
      real*8 root1, root2, root3

!  Locals

      real*8 q, r, m
      real*8 th
      real*8 s,t
      real*8, parameter :: pi = pi_val

      q = (a*a - 3.d0*b)/9.d0
      r = (2.d0*a*a*a - 9.d0*a*b + 27.d0*c)/54.d0
      m = r*r - q*q*q
    
      if (m .lt. 0.d0) then !! three real roots 
        th = acos( r / sqrt(q*q*q) )
        root1 = -(2.d0*sqrt(q)*cos( th         /3.d0))-a/3.d0
        root2 = -(2.d0*sqrt(q)*cos((th+2.d0*pi)/3.d0))-a/3.d0
        root3 = -(2.d0*sqrt(q)*cos((th-2.d0*pi)/3.d0))-a/3.d0
        write(*,*) "three real roots", root1, root2, root3
        stop
      else !! one real root
        if(r .gt. 0.d0) then
          s = -( r+sqrt(m))**(1.d0/3.d0)
        else
          s =  (-r+sqrt(m))**(1.d0/3.d0)
        endif
        t = q/s
        root = s+t-a/3.d0
      endif

      return
      end
